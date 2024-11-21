CC=gcc
CFLAGS=-c -Wall -g
LDFLAGS=-ljpeg -lrt -lpthread -lm
SOURCES= mandel.c jpegrw.c mandelmovie.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLES=mandel mandelmovie

all: $(SOURCES) $(EXECUTABLES)

# pull in dependency info for *existing* .o files
-include $(OBJECTS:.o=.d)

mandel: mandel.o jpegrw.o
	$(CC) mandel.o jpegrw.o $(LDFLAGS) -o $@

mandelmovie: mandelmovie.o jpegrw.o
	$(CC) mandelmovie.o jpegrw.o $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@
	$(CC) -MM $< > $*.d

clean:
	rm -rf $(OBJECTS) $(EXECUTABLES) *.d