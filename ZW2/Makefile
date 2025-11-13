CC =g++
SOURCES = ZW2.c 

ZW2: $(SOURCES)
	$(CC) $(SOURCES) -o $@ -L/usr/local/lib -lsodium -lgmpxx -lgmp -std=c++11 -mcmodel=large -O3
	
.PHONY: clean

clean:
	-rm main

