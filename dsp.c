#include "dsp.h"
#include <nvcom_01.h>

#define DSP_MEM_ADDR(x) ((unsigned int)&x - 0xb8400000)>>2
#define DSP_PTR_ADDR(x) ((unsigned int)x - 0xb8400000)>>2


void init_dsp(int dspNum)
{
	DCSR(dspNum) = 0;
	SR(dspNum) = 0;

	/*A0(dspNum) = DSP_MEM_ADDR(g_nRows1);
	A1(dspNum) = DSP_MEM_ADDR(g_nColumns1);
	A2(dspNum) = DSP_MEM_ADDR(g_nColumns2);
	A0(dspNum) = DSP_PTR_ADDR(g_pReal1);
	A1(dspNum) = DSP_PTR_ADDR(g_pReal2);
	A2(dspNum) = DSP_PTR_ADDR(g_pReal3);
	A3(dspNum) = DSP_MEM_ADDR(g_Factor);*/
}

void run_dsp(int dspNum, unsigned int command)
{
	PC(dspNum) = command;
	DCSR(dspNum) = 0x4000;
	while( !(QSTR_DSP & (1<<3)) ) ;

}
