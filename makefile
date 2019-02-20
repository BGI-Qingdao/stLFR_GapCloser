PROG=GapCloser
all: PROG

PROG:
	@cd Release;make clean;make;cd ..;
	@rm -rf bin;mkdir bin;ln -s ../Release/$(PROG) bin/;

debug:
	@cd Debug;make clean;make;cd ..;
	@rm -rf bin;mkdir bin;ln -s ../Debug/$(PROG) bin/;

clean:
	@cd Release;make clean;cd ..;
	@cd Debug;make clean;cd ..;
	@rm -rf bin;
