/* _qhullmodule.c - Python module for Delaunay triangulations
 *
 * Copyright 2004-2006 Floris Bruynooghe
 *
 * This file is part of Delny.
 *
 * Delny is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Delny is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Delny; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 *
 * Authors: Floris Bruynooghe (flub) <floris.bruynooghe@gmail.com>
 */


#include <Python.h>

#include <stdlib.h>
#include <stdio.h>

#include <numpy/arrayobject.h>
#include <steps/thirdparty/qhull/qhull.h>
#include <steps/thirdparty/qhull/qset.h>		/* for FOREACHneighbor_() */
#include <steps/thirdparty/qhull/poly.h>		/* for qh_vertexneighbors() */


//#define QHULL_31_COMP=1
#ifdef QHULL_31_COMP
char qh_version[] = "_qhullmodule dummy version";
#endif


/* Structures */
struct facet_lists {
	PyObject *facets;
	PyObject *facets_by_index;
};


/* Fucntion declarations */
static PyObject *qhull_delny(PyObject *dummy, PyObject *args);
static void result_qhull_failure(int qhull_exitcode);
static PyObject *result_qhull_success(void);
static struct facet_lists *build_facets(void);
static PyObject *facet_to_pylist(const facetT* facet);
static PyObject *facet_to_pylist_by_index(const facetT* facet);
static PyObject *build_neighbours(void);
static PyObject *vertex_find_neighbours(const vertexT *vertex);
static PyObject *vertex_to_pytuple(const vertexT* vertex);


/* The doc string for the delny function. */
PyDoc_STRVAR(qhull_delny__doc__,
"  delny(nodes) --> (neighbours, facets, facets_by_index)\n\
\n\
A node is represented as a sequence of it's coordinates.  When this\n\
function returns a node this sequence will always be a tuple, however\n\
in it's argument it can be any sequence including Numarray objects.\n\
\n\
points: a sequence type (including Numarray object) of the nodes.\n\
  This generally takes the form of:\n\
  [(x,y,z,...),(x,y,z,...),(x,y,z,...),...]\n\
\n\
neighbours: the neighbours as a dictionary.  Every input point is a\n\
  key and its value is a list of the neighbouring nodes.\n\
\n\
facets: this is a list of facets (triangles in 2D).  Every facet is\n\
  represented by a list of nodes.\n\
\n\
facets_by_index: same as `facets' but uses the index of the position each\n\
  node had in the input list `nodes' to represent a node, not the\n\
  coordinates of the node.");


/** The public Python call.
 *
 * @note See the docstring for documentation.
 */
static PyObject *
qhull_delny(PyObject *dummy, PyObject *args)
{
	PyObject *result = NULL;
	PyObject *input = NULL;	/* input data as we get it from Python */
	PyArrayObject *array = NULL; /* the same input data but converted */
	coordT *points;		/* 1D array of points, dim places a point */
	boolT ismalloc = False;	/* True if qhull should free points in
				   qh_freeqhull() or reallocation */
	int dim;		/* dimension */
	int numpoints;		/* the number of points */
	int exitcode;		/* 0 if no error from qhull */
	char flags[] = "qhull d Qbb Qt"; /* flags for qhull */

	/* Process input. */
	if (!PyArg_ParseTuple (args, "O", &input))
		return NULL;
	array = (PyArrayObject *)
		PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 2, 2);
	if (array == NULL) {
		Py_DECREF(input);
		return NULL;
	}
	dim = array->dimensions[1];
	numpoints = array->dimensions[0];
	points = (coordT *) array->data;

	/* Call qhull */
	exitcode = qh_new_qhull(dim, numpoints, points,
				ismalloc, flags, NULL, stderr);
	if (!exitcode)
		result = result_qhull_success();
	else
		result_qhull_failure(exitcode);

	/* Free the memory */
	qh_freeqhull(qh_ALL);
	Py_DECREF(array);

	/* This is still NULL when qhull failed */
	return result;
} /* qhull_delny() */


/** Set Python exceptions for when qhull failed.
 *
 * @param qhull_exitcode is the exit code returned by the qhull call.
 *
 * @note Maybe there should be specific errors for delny instead of using
 * standard ones?
 */
static void
result_qhull_failure(int qhull_exitcode)
{
	switch (qhull_exitcode) {
	case qh_ERRinput:
		PyErr_BadInternalCall ();
		break;
	case qh_ERRsingular:
		PyErr_SetString(PyExc_ArithmeticError,
				"qhull singular input data");
		break;
	case qh_ERRprec:
		PyErr_SetString(PyExc_ArithmeticError,
				"qhull precision error");
		break;
	case qh_ERRmem:
		PyErr_NoMemory();
		break;
	case qh_ERRqhull:
		PyErr_SetString(PyExc_StandardError,
				"qhull internal error");
		break;
	}
} /* result_qhull_failure() */


/** Create the result of the qhull call.
 *
 * This funcion assumes the qhull call was successful and that all the qhull
 * information can be found under the global `qh' macro/structure.
 *
 * @result New reference to tuple of the neighbours, facets and
 * facets_by_index in that order.
 */
static PyObject *
result_qhull_success(void)
{
	PyObject *result = NULL; /* tuple to return */
	PyObject *neighbours = NULL; /* PyList of neighbours */
	struct facet_lists *facets; /* two PyLists of the facets */

	neighbours = build_neighbours();
	facets = build_facets();
	result = PyTuple_New(3);
	if (neighbours == NULL || facets == NULL || result == NULL)
		goto cleanup;
	PyTuple_SetItem (result, 0, neighbours);
	PyTuple_SetItem (result, 1, facets->facets);
	PyTuple_SetItem (result, 2, facets->facets_by_index);

	return result;

  cleanup:
	Py_XDECREF(result);
	Py_XDECREF(neighbours);
	if (facets != NULL)
		PyMem_Free(facets);
	return NULL;
} /* qhull_result_success() */


/** Build two lists of all the facets; one by coordinates, one by index.
 *
 * A facet is represented by a list of tuples, with each tuple being the
 * coordinates of a node in the facet.  Alternatively as just the index of the
 * node that specifies the position in the input list.
 *
 * @return struct facet_lists.  Note that this needs to be PyMem_Free()ed when
 * it is no longer needed.
 */
static struct facet_lists *
build_facets(void)
{
	struct facet_lists *result; /* return value */
	PyObject *facets = NULL; /* PyList of facets */
	PyObject *facets_by_index = NULL; /* PyList of facets by node index */
	PyObject *tmp = NULL;	/* used to fill up result lists */
	facetT *facet;		/* needed by FORALLfacets */
	int func_retval;	/* to check return values of functions */

	facets = PyList_New(0);
	facets_by_index = PyList_New(0);
	if (facets == NULL || facets_by_index == NULL)
		goto cleanup;
	FORALLfacets {
		if (facet->upperdelaunay)
			continue;
		tmp = facet_to_pylist(facet);
		if (tmp == NULL)
			goto cleanup;
		func_retval = PyList_Append(facets, tmp);
		if (func_retval == -1)
			goto cleanup;
		Py_DECREF(tmp);
		tmp = facet_to_pylist_by_index(facet);
		if (tmp == NULL)
			goto cleanup;
		func_retval = PyList_Append(facets_by_index, tmp);
		if (func_retval == -1)
			goto cleanup;
		Py_DECREF(tmp);
	}
	result = (struct facet_lists *)
		PyMem_Malloc(sizeof(struct facet_lists));
	if (result == NULL) {
		Py_DECREF(facets);
		Py_DECREF(facets_by_index);
		return  (struct facet_lists *) PyErr_NoMemory();
	}
	result->facets = facets;
	result->facets_by_index = facets_by_index;
	return result;

  cleanup:
	Py_XDECREF(facets);
	Py_XDECREF(facets_by_index);
	Py_XDECREF(tmp);
	return NULL;
} /* build_factes() */


/** Convert a facet to it's PyList representation.
 *
 * The PyList representation is a list of tuples where the elements of
 * the tuples are the coordinates of the nodes in the facet.
 *
 * @return new reference to PyList of PyTuple's in facet.
 */
static PyObject *
facet_to_pylist(const facetT* facet)
{
	PyObject *facet_as_list = NULL; /* PyList to return */
	PyObject *node = NULL;	/* PyTuple to fill up facet_as_list */
	vertexT *vertex, **vertexp; /* needed by FOREACHvertex_() */
	int func_retval;	/* to check return values of functions */

	facet_as_list = PyList_New(0);
	if (facet_as_list == NULL)
		return NULL;
	FOREACHvertex_(facet->vertices)	{
		node = vertex_to_pytuple(vertex);
		if (node == NULL)
			goto cleanup;
		func_retval = PyList_Append(facet_as_list, node);
		if (func_retval == -1)
			goto cleanup;
		Py_DECREF(node);
	}

	return facet_as_list;

  cleanup:
	Py_XDECREF(facet_as_list);
	Py_XDECREF(node);
	return NULL;
} /* facet_to_pylist() */


/** Convert a facet to it's PyList (by index) representation.
 *
 * This PyList representation of a facet is a list with only the index
 * numbers (integers) of every node in the facet.
 *
 * @result new reference to PyList of node indices in the facet.
 */
static PyObject *
facet_to_pylist_by_index(const facetT* facet)
{
	PyObject *facet_as_list = NULL; /* PyList to return */
	PyObject *node_index = NULL; /* PyInt to fill up facet_as_list */
	vertexT *vertex, **vertexp; /* needed by FOREACHvertex_() */
	int func_retval;	/* to check return values of functions */

	facet_as_list = PyList_New(0);
	if (facet_as_list == NULL)
		return NULL;
	FOREACHvertex_(facet->vertices)	{
		node_index = PyInt_FromLong(qh_pointid(vertex->point));
		if (node_index == NULL)
			goto cleanup;
		func_retval = PyList_Append(facet_as_list, node_index);
		if (func_retval == -1)
			goto cleanup;
		Py_DECREF(node_index);
	}

	return facet_as_list;

  cleanup:
	Py_XDECREF(facet_as_list);
	Py_XDECREF(node_index);
	return NULL;
} /* facet_to_pylist_by_index() */


/** Build a dictionary of all the neighbouring connections of the convex hull.
 *
 * In the Delaunay triangulation (or convex hull) all vertices are connected
 * to their `nearest neighbours' via `edges' or `ridges'.  This function will
 * build an association between a node and all its neighbours so one can
 * easily find the ridges.
 *
 * @return new reference to a PyDict where the key is a node as a PyTuple of
 *   its coordinates and the value is a PyList of the neighbours as PyTuples of
 *   their coordinates.
 */
static PyObject *
build_neighbours(void)
{
	PyObject *neighbours_dict = NULL; /* return value */
	PyObject *node = NULL;	/* to fill the key of the above dict */
	PyObject *neighbours = NULL; /* to fill the value of the above dict */
	vertexT *vertex;	/* needed by FORALLvertices */

	neighbours_dict = PyDict_New();
	if (neighbours_dict == NULL)
		goto cleanup;
	if (!qh VERTEXneighbors)
		qh_vertexneighbors();
	FORALLvertices {
		node = vertex_to_pytuple(vertex);
		if (node == NULL)
			goto cleanup;
		neighbours = vertex_find_neighbours(vertex);
		if (neighbours == NULL)
			goto cleanup;
		PyDict_SetItem(neighbours_dict, node, neighbours);
		Py_DECREF(node);
		Py_DECREF(neighbours);
	}

	return neighbours_dict;

  cleanup:
	Py_XDECREF(neighbours_dict);
	Py_XDECREF(node);
	Py_XDECREF(neighbours);
	return NULL;
} /* build_neighbours() */


/** Find all the neighbours of a certain vertex.
 *
 * @return new reference to PyList of PyTuples (the neighbours).
 */
static PyObject *
vertex_find_neighbours(const vertexT *myvertex)
{
	PyObject *neighbours = NULL; /* list of neighbours to return */
	PyObject *myvertex_as_pytuple = NULL;
	PyObject *neighbour = NULL;	/* to fill neighbours */
	vertexT *vertex, **vertexp; /* needed by FOREACHvertex_() */
	facetT *neighbor, **neighborp; /* needed by FOREACHneigbor_() */

	neighbours = PyList_New(0);
	myvertex_as_pytuple = vertex_to_pytuple(myvertex);
	if (neighbours == NULL || myvertex_as_pytuple == NULL)
		goto cleanup;
	FOREACHneighbor_(myvertex) {
		/* NOTE: `neighbor' is a facet, not a vertex. */
		if (neighbor->upperdelaunay)
			continue;
		/* FIXME: too much indentation? */
		FOREACHvertex_(neighbor->vertices) {
			neighbour = vertex_to_pytuple(vertex);
			if (neighbour == NULL)
				goto cleanup;
			if (PyObject_Compare(myvertex_as_pytuple, neighbour)
			    && !PySequence_Contains(neighbours, neighbour))
				PyList_Append(neighbours, neighbour);
			Py_DECREF(neighbour);
		}
	}

	return neighbours;

  cleanup:
	Py_XDECREF(neighbours);
	Py_XDECREF(myvertex_as_pytuple);
	Py_XDECREF(neighbour);
	return NULL;
} /* vertex_find_neighbours() */


/** Create a PyTuple from a vertexT.
 *
 * @return a new reference to PyTuple of the coordinates.
 */
static PyObject *
vertex_to_pytuple(const vertexT* vertex)
{
	PyObject *tuple = NULL;
	PyObject *tmp = NULL;
	int i;

	tuple = PyTuple_New(qh hull_dim - 1);
	if (tuple == NULL)
		goto cleanup;
	for (i = 0; i < (qh hull_dim - 1); i++) {
		tmp = PyFloat_FromDouble((double)vertex->point[i]);
		if (tmp == NULL)
			goto cleanup;
		PyTuple_SetItem(tuple, i, tmp);
	}

	return tuple;

  cleanup:
	Py_XDECREF(tuple);
	Py_XDECREF(tmp);
	return NULL;
} /* vertex_to_pytuple() */



/* Tell Python wich methods are available in this module. */
static PyMethodDef QhullMethods[] = {
	{"delny", qhull_delny, METH_VARARGS, qhull_delny__doc__},
	{NULL, NULL, 0, NULL}
};


/* Initialise the module. */
PyMODINIT_FUNC
init_qhull(void)
{
	(void) Py_InitModule("_qhull", QhullMethods);
	import_array()
}



/**** For Emacs
 * Local Variables: 
 * c-file-style: "python"
 * End:
 ****/
