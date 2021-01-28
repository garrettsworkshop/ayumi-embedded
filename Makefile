ayumi.o: ayumi.c ayumi.h
	arm-none-eabi-gcc -mcpu=cortex-m4 -mfloat-abi=hard -Ofast -c ayumi.c

.PHONY: clean
clean:
	rm *.o
