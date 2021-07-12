#ifndef EMP_OT_H__
#define EMP_OT_H__
#include <emp-tool/emp-tool.h>

namespace emp
{

	template <typename T>
	class OT
	{
	public:
		virtual void send(const block *data0, const block *data1, int length) = 0;
		virtual void recv(block *data, const bool *b, int length) = 0;
		virtual ~OT()
		{
		}
	};

/*
*
*special version of OT
*for OT Extension with same choice bit
*/

	template <typename T>
	class OTSameChoiceBit
	{
	public:
		virtual void send() = 0;
		virtual void recv() = 0;
		virtual ~OTSameChoiceBit()
		{
		}
	};

/*
 * 
 * asynchronous ot extension
 */

	template <typename T>
	class OTAsyn
	{
	public:
		virtual void send(const block *data0, const block *data1, int length) = 0;
		virtual void recv(block *data, const bool *b, int length) = 0;
		virtual ~OTAsyn()
		{
		}
	};

//second paper
/*
 * 
 * compress/decompress message in ot extension
 */

	template <typename T>
	class OTWithCom
	{
	public:
		virtual void send(const unsigned int *data0, const unsigned int *data1, int length) = 0;
		virtual void recv(unsigned int *msg, const bool *r, int length) = 0;
		virtual ~OTWithCom()
		{
		}
	};

} // namespace emp
#endif
