/*
 ___license_placeholder___
 */

// Macros to work around improper numpy value conversion support in NumPy 1.4.1

%define UNCHECKED_STL_SEQ_CONVERT(stl_container,insert_method,value_converter)
%typemap(in) const stl_container & {
    PyObject *seq=PySequence_Fast($input, "expected a sequence");
    Py_ssize_t len=PySequence_Size($input);

    stl_container *c=new stl_container();
    for (Py_ssize_t i=0;i<len;++i) 
        c->insert_method(static_cast<stl_container::value_type>(value_converter(PySequence_Fast_GET_ITEM(seq,i))));

    $1=c;
}

%typemap(freearg) const stl_container & {
    if ($1) delete $1;
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const stl_container & {
    $1 = PySequence_Check($input);
}
%enddef

%define UNCHECKED_STL_DICT_CONVERT(stl_container,insert_method,key_converter,value_converter)
%typemap(in) const stl_container & {
    if (!PyDict_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "expected a dictionary");
        return NULL;
    }

    PyObject *keys=PyDict_Keys($input);
    PyObject *values=PyDict_Values($input);
    Py_ssize_t len=PySequence_Size(keys);

    stl_container *c=new stl_container();
    for (Py_ssize_t i=0;i<len;++i) 
        c->insert_method(stl_container::value_type(key_converter(PySequence_Fast_GET_ITEM(keys,i)),value_converter(PySequence_Fast_GET_ITEM(values,i))));

    $1=c;
}

%typemap(freearg) const stl_container & {
    if ($1) delete $1;
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER) const stl_container & {
    $1 = PyDict_Check($input);
}
%enddef

