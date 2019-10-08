#include "dsp.h"
#include "settings.h"

#ifdef DSP_OPTIMIZATION
#include "nvcom_01.h"

#define DSP_MEM_ADDR(x) ((unsigned int)&x - 0xb8400000)>>2
#define DSP_PTR_ADDR(x) ((unsigned int)x - 0xb8400000)>>2

void init_dsp(int dspNum)
{
	DCSR(dspNum) = 0;
	SR(dspNum) = 0;
}

void run_dsp(int dspNum, unsigned int command)
{
	PC(dspNum) = command;
	DCSR(dspNum) = 0x4000;
	while( !(QSTR_DSP & (1<<3)) ) ;
}
#endif // DSP_OPTIMIZATION

#ifndef WIN32
inline unsigned int GetCP0_Count()
{
	unsigned int result;
	asm volatile ("mfc0 %0, $9" :"=r"(result));
	return result;
}
#else // WIN32
	unsigned int GetCP0_Count() { return 0; }
#endif // WIN32
