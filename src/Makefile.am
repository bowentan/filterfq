# CPPFLAGS = -g -std=c++11 -pthread
# LDFLAGS += -Llib
AM_LDFLAGS = -Llib
AM_CPPFLAGS = -g -std=c++11 -static -I../include

bin_PROGRAMS = filterfq
filterfq_SOURCES = filterfq.cpp command_options.cpp fastq_filter.cpp quality_system.cpp
filterfq_LDFLAGS = -L../lib -static
filterfq_LDADD = -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -lboost_thread -lpthread -lz
