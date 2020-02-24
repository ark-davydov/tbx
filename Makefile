
all:
	cd src ; make

clean:	cleanscr cleansrc

cleansrc:
	cd ./src ; make cleanall
cleanscr:
	cd scratch; rm -r ./*

