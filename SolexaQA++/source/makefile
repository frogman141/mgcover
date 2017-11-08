ifeq ($(OS),Windows_NT)
        CC=g++
        CFLAGS=-I.
        BOOST_PATH ?=C:\Boost\include\boost-1_59
        ZLIB_PATH ?= C:\MinGW
        DEPS= kseq.h R_codes.h
        CXXFLAGS=-std=c++1y
        LDFLAGS=-lm
        LIBS= -lboost_filesystem-mgw61-mt-1_59 -lboost_system-mgw61-mt-1_59 -lboost_regex-mgw61-mt-1_59 -lboost_iostreams-mgw61-mt-1_59 -lz
        %.o: %.c $(DEPS)
	        $(CC) -c -o $@ $< $(CFLAGS)

        SolexaQA++: main_win.cpp
	        $(CC) -static main_win.cpp -o SolexaQA++ $(CXXFLAGS) -I. -I$(BOOST_PATH) -I$(ZLIB_PATH)/include -L$(ZLIB_PATH)/lib -L$(BOOST_PATH)/../../lib $(LDFLAGS) $(LIBS)

else
        UNAME_S := $(shell uname -s)
        ifeq ($(UNAME_S),Linux)
                CC=g++
                CFLAGS=-I.
                BOOST_PATH ?= /usr/
                DEPS= kseq.h R_codes.h
        	CXXFLAGS=-std=c++11
                LDFLAGS=-lm
                LIBS=-lz -lboost_filesystem -lboost_system -lboost_regex -lboost_iostreams -lpthread
                %.o: %.c $(DEPS)
	                $(CC) -c -o $@ $< $(CFLAGS)

                SolexaQA++: main.cpp 
	            	$(CC) -static main.cpp -o SolexaQA++ -fstack-protector-all -I$(BOOST_PATH)/include $(CFLAGS) -L$(BOOST_PATH)/lib $(LDFLAGS) $(CXXFLAGS) $(LIBS)
        endif
        ifeq ($(UNAME_S),Darwin)
                CC=g++
                DEPS= kseq.h R_codes.h
                BOOST_PATH ?=/usr/local/
                CXXFLAGS= -std=c++11 -stdlib=libc++ -mmacosx-version-min=10.11
                LDFLAGS=-lm
                LIBS=-lz -lpthread
                
                ifeq ($(LIBBOOST_FS),)
                        ifneq ($(wildcard $(BOOST_PATH)lib/libboost_filesystem.a),)
                            LIBBOOST_FS = $(BOOST_PATH)lib/libboost_filesystem.a
                        endif
                else 
                        LIBBOOST_FS = $(BOOST_PATH)lib/libboost_filesystem-mt.a
                endif
                
                ifeq ($(LIBBOOST_SYS),)
                        ifneq ($(wildcard $(BOOST_PATH)lib/libboost_system.a),)
                            LIBBOOST_SYS = $(BOOST_PATH)lib/libboost_system.a
                        endif
                else 
                        LIBBOOST_SYS = $(BOOST_PATH)lib/libboost_system-mt.a
                endif
                
                ifeq ($(LIBBOOST_REGEX),)
                        ifneq ($(wildcard $(BOOST_PATH)lib/libboost_regex.a),)
                            LIBBOOST_REGEX = $(BOOST_PATH)lib/libboost_regex.a
                        endif
                else 
                        LIBBOOST_REGEX = $(BOOST_PATH)lib/libboost_regex-mt.a
                endif
                                
                ifeq ($(LIBBOOST_IO),)
                        ifneq ($(wildcard $(BOOST_PATH)lib/libboost_iostreams.a),)
                            LIBBOOST_IO = $(BOOST_PATH)lib/libboost_iostreams.a
                        endif
                else 
                        LIBBOOST_IO = $(BOOST_PATH)lib/libboost_iostreams-mt.a
                endif
                                
                SolexaQA++: main.cpp
			$(CC) main.cpp -o SolexaQA++ $(LIBBOOST_FS) $(LIBBOOST_SYS) $(LIBBOOST_REGEX) $(LIBBOOST_IO) -I. -I/usr/local/include/boost -I/usr/local/include/ -L/usr/local/lib $(LDFLAGS) $(CXXFLAGS) $(LIBS)

endif
endif
