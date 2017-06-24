cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream& write(const char*, int) except +
    cdef cppclass istream:
        istream& read(const char*, int) except +

# obviously std::ios_base isn't a namespace, but this lets
# Cython generate the connect C++ code
cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef cppclass streambuf:
        pass
    cdef cppclass iostream(istream):
        iostream(streambuf*)
    cdef open_mode binary

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        # constructors
        ofstream(const char*) except +
        ofstream(const char*, open_mode) except+
    
    cdef cppclass ifstream(istream):
        # constructors
        ifstream(const char*) except +
        ifstream(const char*, open_mode) except+
    
    cdef cppclass filebuf:
        pass
        
    cdef cppclass fstream():
        fstream()
        fstream(const char*, openmode)
        fstream(const char&, openmode)
        void close()
        bool is_open()
        void open(const char*, openmode)
        void open(const char&, openmode)
        filebuf* rdbuf() const
        
        