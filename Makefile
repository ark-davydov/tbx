
all:
	cd src ; make

clean:	cleansrc cleanscr

cleansrc:
	cd ./src ; make cleanall
cleanscr:
	cd scratch; rm -r ./*

