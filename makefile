
CC		= g++
SRC_DIR		= src
EXEC_DIR	= bin
SRC		= main_kcmbt.cc
EXEC		= kcmbt
DUMP_SRC	= decode.cc
DUMP_EXEC	= kcmbt_dump
CFLAGS		= -Wno-unused-result -std=c++11 -O3
CLINK		= -lz 

all: kcmbt kcmbt_dump

kcmbt:
	mkdir -p bin
	$(CC) $(CFLAGS) $(SRC_DIR)/$(SRC) -o $(EXEC_DIR)/$(EXEC) $(CLINK)
	
kcmbt_dump:
	$(CC) $(CFLAGS) $(SRC_DIR)/$(DUMP_SRC) -o $(EXEC_DIR)/$(DUMP_EXEC) $(CLINK)
	
	
clean:
	rm -f $(EXEC_DIR)/$(EXEC) $(EXEC_DIR)/$(DUMP_EXEC)


