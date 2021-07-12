#ifndef EMP_COT_H__
#define EMP_COT_H__
#include "emp-ot/ot.h"

namespace emp
{
	const static int ot_bsize = 8;
	template <typename T>
	class COT : public OT<T>
	{
	public:
		T *io = nullptr;
		MITCCRH<ot_bsize> mitccrh;
		block Delta;
		virtual void send_cot(block *data0, int length) = 0;
		virtual void recv_cot(block *data, const bool *b, int length) = 0;
		void send(const block *data0, const block *data1, int length) override
		{
			block *data = new block[length];
			send_cot(data, length);
			mitccrh.setS(zero_block);
			io->flush();
			block pad[2 * ot_bsize];

			auto start = clock_start();
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data[j];
					pad[2 * (j - i) + 1] = data[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);

				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = pad[2 * (j - i)] ^ data0[j];
					pad[2 * (j - i) + 1] = pad[2 * (j - i) + 1] ^ data1[j];
				}
				io->send_data(pad, 2 * sizeof(block) * min(ot_bsize, length - i));
			}
			//std::cout << "aes time cost " << time_from(start) * 1.0 / 1000000 << endl;

			delete[] data;
		}

		void recv(block *data, const bool *r, int length) override
		{
			recv_cot(data, r, length);
			mitccrh.setS(zero_block);
			io->flush();

			block res[2 * ot_bsize];
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				io->recv_data(res, 2 * sizeof(block) * min(ot_bsize, length - i));
				for (int j = 0; j < ot_bsize and j < length - i; ++j)
				{
					data[i + j] = res[2 * j + r[i + j]] ^ pad[j];
				}
			}
		}

		void send_rot(block *data0, block *data1, int length)
		{
			send_cot(data0, length);
			io->flush();
			block pad[ot_bsize * 2];
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data0[j];
					pad[2 * (j - i) + 1] = data0[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					data0[j] = pad[2 * (j - i)];
					data1[j] = pad[2 * (j - i) + 1];
				}
			}
		}

		void recv_rot(block *data, const bool *r, int length)
		{
			recv_cot(data, r, length);
			io->flush();
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				memcpy(data + i, pad, min(ot_bsize, length - i) * sizeof(block));
			}
		}
	};

	/*
*
*special version of COT
*for OT Extension with same choice bit
*/
	template <typename T>
	class COTSameChoiceBit : public OTSameChoiceBit<T>
	{
	public:
		T *io = nullptr;
		MITCCRH<ot_bsize> mitccrh;
		block Delta;
		block *otex_data;
		virtual void send_cot(block *data0, int length) = 0;
		virtual void recv_cot(block *data, const bool *b, int length) = 0;

		void set_up_sender(long lenght)
		{
			otex_data = new block[lenght];
			send_cot(otex_data, lenght);
			io->flush();
		}

		void set_up_receiver(const bool *r, long lenght)
		{
			otex_data = new block[lenght];
			recv_cot(otex_data, r, lenght);
			io->flush();
		}

		void encrypt(block *m1_block, block *m0_block, uint8_t *ot_pad_shift, long m1_lenght, long m0_lenght)
		{
			mitccrh.setS(zero_block);

			block pad_m[ot_bsize];
			long pad_index = -1;
			//encrypt m1
			for (long i = 0; i < m1_lenght;)
			{
				for (int j = 0; j < ot_bsize; j++)
				{
					pad_index += ot_pad_shift[0];
					ot_pad_shift += 1;
					pad_m[j] = otex_data[pad_index] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 1>(pad_m);

				for (int j = 0; j < ot_bsize; j++)
					m1_block[j] = pad_m[j] ^ m1_block[j];
				i += ot_bsize;
				m1_block += ot_bsize;
			}

			//encrypt m0
			for (long i = 0; i < m0_lenght;)
			{
				for (int j = 0; j < ot_bsize; j++)
					pad_m[j] = otex_data[i + j];
				mitccrh.hash<ot_bsize, 1>(pad_m);

				for (int j = 0; j < ot_bsize; j++)
					m0_block[j] = pad_m[j] ^ m0_block[j];
				i += ot_bsize;
				m0_block += ot_bsize;
			}

			delete[] otex_data;
		}

		void decrypt(block *m1_block, block *m0_block, uint8_t *ot_pad_shift, long m1_lenght, long m0_lenght)
		{
			mitccrh.setS(zero_block);

			block pad_m[ot_bsize];
			//decrypt m1
			long pad_index = -1;
			for (long i = 0; i < m1_lenght;)
			{
				for (int j = 0; j < ot_bsize; j++)
				{
					pad_index += ot_pad_shift[0];
					ot_pad_shift += 1;
					pad_m[j] = otex_data[pad_index];
				}
				mitccrh.hash<ot_bsize, 1>(pad_m);

				for (int j = 0; j < ot_bsize; j++)
					m1_block[j] = m1_block[j] ^ pad_m[j];
				m1_block += ot_bsize;
				i += ot_bsize;
			}

			//decrypt m0
			for (long i = 0; i < m0_lenght;)
			{
				for (int j = 0; j < ot_bsize; j++)
					pad_m[j] = otex_data[i + j];
				mitccrh.hash<ot_bsize, 1>(pad_m);
				for (int j = 0; j < ot_bsize; j++)
					m0_block[j] = m0_block[j] ^ pad_m[j];

				m0_block += ot_bsize;
				i += ot_bsize;
			}

			delete[] otex_data;
		}

		void send_encrypted_data(block *m1_block, block *m0_block, long m1_lenght, long m0_lenght)
		{
			//send 1GB each io round  67108864 * 16 = 1GB
			int block_send_size = 67108864;
			for (long i = 0; i < m1_lenght;)
			{
				block_send_size = std::min((long)67108864, m1_lenght - i);
				io->send_data(m1_block + i, sizeof(block) * block_send_size);
				i += block_send_size;
			}
			for (long i = 0; i < m0_lenght;)
			{
				block_send_size = std::min((long)67108864, m0_lenght - i);
				io->send_data(m0_block + i, sizeof(block) * block_send_size);
				i += block_send_size;
			}
			io->flush();
		}

		void recv_encrypted_data(block *m1_block, block *m0_block, long m1_lenght, long m0_lenght)
		{
			//recv 1GB each io round  67108864 * 16 = 1GB
			int block_recv_size = 67108864;

			for (long i = 0; i < m1_lenght;)
			{
				block_recv_size = std::min((long)67108864, m1_lenght - i);
				io->recv_data(m1_block + i, sizeof(block) * block_recv_size);
				i += block_recv_size;
			}
			for (long i = 0; i < m0_lenght;)
			{
				block_recv_size = std::min((long)67108864, m0_lenght - i);
				io->recv_data(m0_block + i, sizeof(block) * block_recv_size);
				i += block_recv_size;
			}
			io->flush();
		}

		void send() override
		{
		}

		void recv() override
		{
		}

		void send_rot(block *data0, block *data1, int length)
		{
			send_cot(data0, length);
			io->flush();
			block pad[ot_bsize * 2];
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data0[j];
					pad[2 * (j - i) + 1] = data0[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					data0[j] = pad[2 * (j - i)];
					data1[j] = pad[2 * (j - i) + 1];
				}
			}
		}

		void recv_rot(block *data, const bool *r, int length)
		{
			recv_cot(data, r, length);
			io->flush();
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				memcpy(data + i, pad, min(ot_bsize, length - i) * sizeof(block));
			}
		}
	};

	/*
*
*asynchronous OT extension
*/
	template <typename T>
	class COTAsyn : public OTAsyn<T>
	{
	public:
		T *io = nullptr;
		MITCCRH<ot_bsize> mitccrh;
		block Delta;
		block *pad_data;
		long all_length, index;
		bool *cot_r;
		virtual void send_cot(block *data0, int length) = 0;
		virtual void recv_cot(block *data, const bool *b, int length) = 0;

		void asyn_setup_send(long all_length)
		{
			this->all_length = all_length;
			pad_data = new block[all_length];
			send_cot(pad_data, all_length);
			io->flush();
			index = 0;
			std::cout << "asyn_setup_send ok" << endl;
		}

		void asyn_setup_recv(long all_length)
		{
			cot_r = new bool[all_length];
			for (long i = 0; i < all_length; i++)
				cot_r[i] = rand() % 2;

			this->all_length = all_length;
			pad_data = new block[all_length];
			recv_cot(pad_data, cot_r, all_length);
			io->flush();
			index = 0;
			std::cout << "asyn_setup_recv ok" << endl;
		}

		void print128_num(__m128i *var, int lenght)
		{
			for (int i = 0; i < lenght; i++)
			{
				uint16_t *ptr = (uint16_t *)&(var[i]);
				for (int j = 0; j < sizeof(block) / sizeof(uint16_t); j++)
					printf("%4x ", ptr[j]);
				printf("\n");
			}
		}

		void send(const block *data0, const block *data1, int length) override
		{
			if (index + length > all_length)
			{
				std::cout << "re setup" << endl;
				delete[] pad_data;
				asyn_setup_send(1000 * 10000);
			}

			bool *rd = new bool[length];
			io->recv_data(rd, length);

			block *data = pad_data + index;
			mitccrh.setS(zero_block);
			block pad[2 * ot_bsize];

			auto start = clock_start();
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data[j];
					pad[2 * (j - i) + 1] = data[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);

				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					if (rd[j])
					{
						block tmp = pad[2 * (j - i) + 1] ^ data0[j];
						pad[2 * (j - i) + 1] = pad[2 * (j - i)] ^ data1[j];
						pad[2 * (j - i)] = tmp;
					}
					else
					{
						pad[2 * (j - i)] = pad[2 * (j - i)] ^ data0[j];
						pad[2 * (j - i) + 1] = pad[2 * (j - i) + 1] ^ data1[j];
					}
				}
				io->send_data(pad, 2 * sizeof(block) * min(ot_bsize, length - i));
			}

			delete[] rd;
			index += length;
		}

		void recv(block *data, const bool *r, int length) override
		{
			if (index + length > all_length)
			{
				std::cout << "re setup" << endl;
				delete[] pad_data;
				delete[] cot_r;
				asyn_setup_recv(1000 * 10000);
			}
			for (long i = 0; i < length; i++)
				cot_r[i + index] = cot_r[i + index] ^ r[i];
			io->send_data(cot_r + index, length);

			mitccrh.setS(zero_block);

			block res[2 * ot_bsize];
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, pad_data + index + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				io->recv_data(res, 2 * sizeof(block) * min(ot_bsize, length - i));
				for (int j = 0; j < ot_bsize and j < length - i; ++j)
				{
					data[i + j] = res[2 * j + r[i + j]] ^ pad[j];
				}
			}

			index += length;
		}

		void send_rot(block *data0, block *data1, int length)
		{
			send_cot(data0, length);
			io->flush();
			block pad[ot_bsize * 2];
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data0[j];
					pad[2 * (j - i) + 1] = data0[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					data0[j] = pad[2 * (j - i)];
					data1[j] = pad[2 * (j - i) + 1];
				}
			}
		}

		void recv_rot(block *data, const bool *r, int length)
		{
			recv_cot(data, r, length);
			io->flush();
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				memcpy(data + i, pad, min(ot_bsize, length - i) * sizeof(block));
			}
		}
	};

	//second paper
	/*
*
*compress/decompress message in ot extension
*/
	template <typename T>
	class COTWithCom : public OTWithCom<T>
	{
	public:
		T *io = nullptr;
		MITCCRH<ot_bsize> mitccrh;
		block Delta;
		virtual void send_cot(block *data0, int length) = 0;
		virtual void recv_cot(block *data, const bool *b, int length) = 0;

		int getComLen(int length)
		{
			int numel = length / 32;
			return numel * 17;
		}

		void otCompress(unsigned int *comed, unsigned int *data, int length)
		{
			int comIndex = 0;

			for (int i = 0; i < length / 32; i++)
			{
				for (int j = 0; j <= 16; j++)
				{
					if (j == 0 || j == 16)
					{
						comed[comIndex] = data[j];
					}
					else
					{
						unsigned int t0 = (data[j] >> j) << j;
						unsigned int t1 = data[32 - j] >> (32 - j);
						comed[comIndex] = t0 | t1;
					}

					comIndex++;
				}
				data += 32;
			}
		}

		void otDecompress(unsigned int *decomed, unsigned int *data, int length)
		{
			int decomIndex = 0;
			int itemLen = (2 + (32 - 2) / 2);

			for (int i = 0; i < length / 32; i++)
			{
				for (int j = 0; j <= 16; j++)
				{
					if (j == 0 || j == 16)
					{
						decomed[decomIndex + j] = data[j];
					}
					else
					{
						int t0 = data[j];
						int t1 = t0 << (32 - j);
						decomed[decomIndex + j] = (t0 >> j) << j;
						decomed[decomIndex + 32 - j] = t1;
					}
				}
				decomIndex += 32;
				data += itemLen;
			}
		}

		void send(const unsigned int *data0, const unsigned int *data1, int length) override
		{
			block *data = new block[length];
			send_cot(data, length);
			mitccrh.setS(zero_block);
			io->flush();

			block pad[2 * ot_bsize];
			block *mask = new block[length * 2];
			//auto start = clock_start();
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data[j];
					pad[2 * (j - i) + 1] = data[j] ^ Delta;
				}

				mitccrh.hash<ot_bsize, 2>(pad);
				memcpy(mask + i * 2, pad, 2 * sizeof(block) * min(ot_bsize, length - i));
			}
			//std::cout << "aes time cost " << time_from(start) * 1.0 / 1000000 << endl;
			delete[] data;

			unsigned int *maskUint = (unsigned int *)mask;
			unsigned int *m0Masked = new unsigned int[length];
			unsigned int *m1Masked = new unsigned int[length];
			for (int i = 0; i < length; i++)
			{
				m0Masked[i] = maskUint[8 * i] ^ data0[i];
				m1Masked[i] = maskUint[8 * i + 4] ^ data1[i];
			}

			//compress
			int comlen = getComLen(length);
			unsigned int *m0Comed = new unsigned int[comlen];
			unsigned int *m1Comed = new unsigned int[comlen];
			otCompress(m0Comed, m0Masked, length);
			otCompress(m1Comed, m1Masked, length);
            
			//send ot msg
			io->send_data(m0Comed, sizeof(int) * comlen);
			io->send_data(m1Comed, sizeof(int) * comlen);
			io->flush();

			delete[] mask;
			delete[] m0Masked;
			delete[] m1Masked;
			delete[] m0Comed;
			delete[] m1Comed;
		}

		void recv(unsigned int *msg, const bool *r, int length) override
		{
			block *data = new block[length];
			recv_cot(data, r, length);
			mitccrh.setS(zero_block);
			io->flush();

			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				memcpy(data + i, pad, min(ot_bsize, length - i) * sizeof(block));
			}
            
			//recv ot msg
			int comlen = getComLen(length);
			unsigned int *maskUint = (unsigned int *)data;
			unsigned int *m0Recved = new unsigned int[comlen];
			unsigned int *m1Recved = new unsigned int[comlen];
			io->recv_data(m0Recved, sizeof(int) * comlen);
			io->recv_data(m1Recved, sizeof(int) * comlen);
			io->flush();
            
			//decompress
			unsigned int *m0Decomed = new unsigned int[length];
			unsigned int *m1Decomed = new unsigned int[length];
			otDecompress(m0Decomed, m0Recved, length);
			otDecompress(m1Decomed, m1Recved, length);

			for (int i = 0; i < length; i++)
			{
				msg[i] = (1 - r[i]) * (m0Decomed[i] ^ maskUint[i * 4]) + r[i] * (m1Decomed[i] ^ maskUint[i * 4]);
				int shift = i % 32;
				msg[i] = (msg[i] >> shift) << (shift);
			}

			delete[] data;
			delete[] m0Recved;
			delete[] m1Recved;
			delete[] m0Decomed;
			delete[] m1Decomed;
		}

		void send_rot(block *data0, block *data1, int length)
		{
			send_cot(data0, length);
			io->flush();
			block pad[ot_bsize * 2];
			for (int i = 0; i < length; i += ot_bsize)
			{
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					pad[2 * (j - i)] = data0[j];
					pad[2 * (j - i) + 1] = data0[j] ^ Delta;
				}
				mitccrh.hash<ot_bsize, 2>(pad);
				for (int j = i; j < min(i + ot_bsize, length); ++j)
				{
					data0[j] = pad[2 * (j - i)];
					data1[j] = pad[2 * (j - i) + 1];
				}
			}
		}

		void recv_rot(block *data, const bool *r, int length)
		{
			recv_cot(data, r, length);
			io->flush();
			block pad[ot_bsize];
			for (int i = 0; i < length; i += ot_bsize)
			{
				memcpy(pad, data + i, min(ot_bsize, length - i) * sizeof(block));
				mitccrh.hash<ot_bsize, 1>(pad);
				memcpy(data + i, pad, min(ot_bsize, length - i) * sizeof(block));
			}
		}
	};

} // namespace emp
#endif
