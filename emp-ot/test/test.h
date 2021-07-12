#include <emp-tool/emp-tool.h>
#include "emp-ot/emp-ot.h"
#include <iostream>
using namespace emp;
using namespace std;


template <typename T>
double test_ot(T * ot, NetIO *io, int party, int length) {
	block *b0 = new block[length], *b1 = new block[length],
	*r = new block[length];
	PRG prg(fix_key);
	prg.random_block(b0, length);
	prg.random_block(b1, length);
	bool *b = new bool[length];
	PRG prg2;
	prg2.random_bool(b, length);

	auto start = clock_start();
	if (party == ALICE) {
		ot->send(b0, b1, length);
	} else {
		ot->recv(r, b, length);
	}
	io->flush();
	long long t = time_from(start);
	std::cout << "time cost " << t * 1.0 / 1000000.0 << endl;

	std::cout << "Tests passed.\t";
	delete[] b0;
	delete[] b1;
	delete[] r;
	delete[] b;
	return t;
}

void print128_num(__m128i var)
{
    uint16_t val[8];
    memcpy(val, &var, sizeof(val));
    printf("%4x %4x %4x %4x %4x %4x %4x %4x \n", 
           val[0], val[1], val[2], val[3], val[4], val[5], 
           val[6], val[7]);
}

void print128_num(__m128i *var ,int lenght)
{
	for (int i = 0; i < lenght; i++)
	{
		uint16_t *ptr = (uint16_t *)&(var[i]);
		for (int j = 0; j < sizeof(block) / sizeof(uint16_t); j++)
			printf("%4x ", ptr[j]);
		printf("\n");
	}
}


template <typename T>
void test_ot_setup(T *ot, NetIO *io, int party, long length)
{
	length = 128;
	block *m0 = new block[length], *m1 = new block[length];
	PRG prg(fix_key);
	prg.random_block(m0, length);
	prg.random_block(m1, length);

	if(party == ALICE)
		ot->asyn_setup_send(1000 * 10000);
	else
		ot->asyn_setup_recv(1000 * 10000);

	std::cout << "m0" << endl;
	//print128_num(m0, length);
	std::cout << "m1" << endl;
	//print128_num(m1, length);
    
	bool *b = new bool[length];
	for (int i = 0; i < length; i++)
		b[i] = rand() % 2;

	auto start = clock_start();

	for (int i = 0; i < 10; i++)
	{
		if (party == ALICE)
		{
			ot->send(m0, m1, length);
		}
		else
		{
			ot->recv(m0, b, length);
		}
		io->flush();

		//std::cout << "recved" << endl;
		//print128_num(m0, 2);
	}

	long long t = time_from(start);
	cout << "time cost " << t * 1.0 / 1000000.0 << endl;

	delete[] m0;
	delete[] m1;
}

template <typename T>
double test_ot_samebit(T *ot, NetIO *io, int party, int length)
{
	int batchblock = 2, itemnum = 4;
	length = batchblock * itemnum;

	block *b0 = new block[itemnum], *b1 = new block[length];
	block *m0 = new block[itemnum], *m1 = new block[length];
	PRG prg(fix_key);
	prg.random_block(b0, itemnum);
	prg.random_block(b1, length);
	bool *b = new bool[itemnum];
	PRG prg2;
	prg2.random_bool(b, itemnum);
    
	//io->counter = 0;
	for (int i = 0; i < 1; i++)
	{
		cout << i << endl;
		auto start = clock_start();

		if (party == ALICE)
		{
			ot->send(b0, b1, length, batchblock);
		}
		else
		{
			ot->recv(m0, m1, b, length, batchblock);
		}
		io->flush();

		//cout << "counter " << io->counter << endl;
		long long t = time_from(start);
		cout << "time cost " << t * 1.0 / 1000000.0 << endl;
	}

	std::cout << "Tests passed.\n";
	delete[] b0;
	delete[] b1;
	delete[] m0;
	delete[] m1;
	return 0;
}

template <typename T>
double test_cot(T *ot, NetIO *io, int party, int length)
{
	block *b0 = new block[length], *r = new block[length];
	bool *b = new bool[length];
	block delta;
	PRG prg;
	prg.random_block(&delta, 1);
	prg.random_bool(b, length);

	io->sync();
	auto start = clock_start();
	if (party == ALICE)
	{
		ot->send_cot(b0, length);
		delta = ot->Delta;
	}
	else
	{
		ot->recv_cot(r, b, length);
	}
	io->flush();
	long long t = time_from(start);
	if (party == ALICE)
	{
		io->send_block(&delta, 1);
		io->send_block(b0, length);
	}
	else if (party == BOB)
	{
		io->recv_block(&delta, 1);
		io->recv_block(b0, length);
		for (int i = 0; i < length; ++i)
		{
			block b1 = b0[i] ^ delta;
			if (b[i])
			{
				if (!cmpBlock(&r[i], &b1, 1))
					error("COT failed!");
			}
			else
			{
				if (!cmpBlock(&r[i], &b0[i], 1))
					error("COT failed!");
			}
		}
	}
	std::cout << "Tests passed.\t";
	io->flush();
	delete[] b0;
	delete[] r;
	delete[] b;
	return t;
}

template <typename T>
double test_rot(T *ot, NetIO *io, int party, int length)
{
	block *b0 = new block[length], *r = new block[length];
	block *b1 = new block[length];
	bool *b = new bool[length];
	PRG prg;
	prg.random_bool(b, length);

	io->sync();
	auto start = clock_start();
	if (party == ALICE)
	{
		ot->send_rot(b0, b1, length);
	}
	else
	{
		ot->recv_rot(r, b, length);
	}
	io->flush();
	long long t = time_from(start);
	if (party == ALICE)
	{
		io->send_block(b0, length);
		io->send_block(b1, length);
	}
	else if (party == BOB)
	{
		io->recv_block(b0, length);
		io->recv_block(b1, length);

		/*
		int64_t *v64val = (int64_t *)&(b0[0]);
		int b_value = v64val[0];
		int v = b_value;
		v64val = (int64_t *)&(b1[0]);
		b_value = b_value xor v64val[0];
		std::cout << endl;
		std::cout << "b" << endl
			 << b_value << endl;
		std::cout << "v" << endl
			 << v << endl;

		v64val = (int64_t *)&(r[0]);
		int u = v64val[0];
		std::cout << "a" << endl
			 << b[0] << endl;
		std::cout << "u" << endl
			 << u << endl;
		std::cout << "u xor v" << endl
				  << (u xor v) << endl;
		*/

		for (int i = 0; i < length; ++i)
		{
			if (b[i])
				assert(cmpBlock(&r[i], &b1[i], 1));
			else
				assert(cmpBlock(&r[i], &b0[i], 1));
		}
	}
	std::cout << "Tests passed.\t";
	io->flush();
	delete[] b0;
	delete[] b1;
	delete[] r;
	delete[] b;
	return t;
}

template <typename T>
double test_rcot(T *ot, NetIO *io, int party, int length, bool inplace)
{
	block *b = nullptr;
	PRG prg;

	io->sync();
	auto start = clock_start();
	uint64_t mem_size;
	if (!inplace)
	{
		mem_size = length;
		b = new block[length];

		// The RCOTs will be generated in the internal buffer
		// then be copied to the user buffer
		ot->rcot(b, length);
	}
	else
	{
		// Call byte_memory_need_inplace() to get the buffer size needed
		mem_size = ot->byte_memory_need_inplace((uint64_t)length);
		b = new block[mem_size];

		// The RCOTs will be generated directly to this buffer
		ot->rcot_inplace(b, mem_size);
	}
	long long t = time_from(start);
	io->sync();
	if (party == ALICE)
	{
		io->send_block(&ot->Delta, 1);
		io->send_block(b, mem_size);
	}
	else if (party == BOB)
	{
		block ch[2];
		ch[0] = zero_block;
		block *b0 = new block[mem_size];
		io->recv_block(ch + 1, 1);
		io->recv_block(b0, mem_size);
		for (size_t i = 0; i < mem_size; ++i)
		{
			b[i] = b[i] ^ ch[getLSB(b[i])];
		}
		if (!cmpBlock(b, b0, mem_size))
			error("RCOT failed");
		delete[] b0;
	}
	std::cout << "Tests passed.\t";
	delete[] b;
	return t;
}
