
all:
	mkdir -p scratch; cp -r examples/* scratch; cd src ; make

clean:	cleanscr cleansrc

cleansrc:
	cd src ; make cleanall

# mkdir is needed, otherwise this command may delete evergything
cleanscr:
	mkdir -p scratch; cd scratch; rm -r ./*

