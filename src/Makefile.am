# CPPFLAGS = -g -std=c++11 -pthread
# LDFLAGS += -Llib
AM_LDFLAGS = -Llib
AM_CPPFLAGS = -g -std=c++11 -pthread -I../include

bin_PROGRAMS = filterfq
filterfq_SOURCES = main.cpp command_options.cpp fastq_filter.cpp quality_system.cpp ../include/command_options.hpp ../include/fastq_filter.hpp ../include/quality_system.hpp
filterfq_LDADD = -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -lboost_thread
