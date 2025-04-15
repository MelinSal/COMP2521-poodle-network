# COMP2521 - Assignment 2

# List all your supporting .c files here. Do NOT include .h files in this list.
# Example: SUPPORTING_FILES = hello.c world.c

SUPPORTING_FILES = 

########################################################################
# !!! DO NOT MODIFY ANYTHING BELOW THIS LINE !!!

CC = clang
CFLAGS = -Wall -Wvla -Werror -gdwarf-4

########################################################################

.PHONY: asan msan nosan

asan: CFLAGS += -fsanitize=address,leak,undefined
asan: all

msan: CFLAGS += -fsanitize=memory,undefined -fsanitize-memory-track-origins
msan: all

nosan: all

########################################################################

.PHONY: all
all: testPoodle

testPoodle: testPoodle.c poodle.c $(SUPPORTING_FILES)
	$(CC) $(CFLAGS) -o testPoodle testPoodle.c poodle.c $(SUPPORTING_FILES)

.PHONY: clean
clean:
	rm -f testPoodle
