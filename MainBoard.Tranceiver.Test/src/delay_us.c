/*
 * delay_us.c
 *
 *  Created on: Nov 27, 2011
 *      Author: Paul Shi
 */
#include "delay_us.h"


void delay_us(unsigned long delay,unsigned long fcpu_hz)
{
  cpu_delay_us(delay, fcpu_hz);
}