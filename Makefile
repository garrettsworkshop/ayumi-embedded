ayumi.o: ayumi.c ayumi.h
	arm-none-eabi-gcc -mcpu=cortex-m4 -mfloat-abi=hard -O3 -c ayumi.c
