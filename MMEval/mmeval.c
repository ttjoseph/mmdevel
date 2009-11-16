/**
 * Python extension to evaluate molecular mechanics energies in a simple
 * manner. That is, with no cutoff, pairlists, or such. These methods
 * take an AmberSystem object.
 *
 * Tom Joseph <thomas.joseph@mssm.edu>
 */
#include <Python.h>

/**
 * Returns a PyObject representing the AmberSystem "blocks"
 * associative array. Used to help in extracting information
 * from Python types that we receive from the interpreter.
 */
PyObject *get_block(PyObject *mol, const char *block_name)
{
  // Extract the blocks object from the AmberSystem
  PyObject *blocks;
  if((blocks = PyObject_GetAttrString(mol, "blocks")) == NULL)
  {
    fprintf(stderr, "get_block: Couldn't find blocks attribute\n");
    return NULL;
  }
  
  // Get the block of the specified name
  PyObject *block;
  if((block = PyDict_GetItemString(blocks, block_name)) == NULL)
  {
    fprintf(stderr, "get_block: Couldn't find block %s\n", block_name);
    return NULL;
  }

  return block;
}

/**
 * Gets an array of floats from an AmberSystem object.
 * Intended to be used for x, y, z coordinate arrays.
 *
 * Equivalent to: mol.attr_name
 */
float *get_float_array(PyObject *mol, const char *attr_name, int *size)
{
  PyObject *data_obj;
  if((data_obj = PyObject_GetAttrString(mol, attr_name)) == NULL)
  {
    fprintf(stderr, "get_float_array: Couldn't find attribute %s\n", attr_name);
    return NULL;
  }
    
  // Extract the list of doubles
  Py_ssize_t length = PySequence_Size(data_obj);
  *size = (int) length;
  float *data = (float*) malloc(sizeof(float) * (int) length);
  for(Py_ssize_t i = 0; i < length; i++)
    data[i] = (float) PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(data_obj, i));
  
  return data;  
}

/**
 * Returns an AmberSystem block as an array of ints.
 *
 * You have to free the thing you get back from this method!
 * Equivalent to: mol.blocks[block_name]
 */
int *get_block_as_int_array(PyObject *mol, const char *block_name, int *size)
{
  PyObject *block;
  if((block = get_block(mol, block_name)) == NULL)
  {
    // fprintf(stderr, "get_block_as_int_array: Couldn't find block %s\n", block_name);
    return NULL;
  }
  // Extract the list of integers from the AmberSystem
  Py_ssize_t block_length = PySequence_Size(block);
  *size = (int) block_length;
  int *data = (int*) malloc(sizeof(int) * (int) block_length);
  for(Py_ssize_t i = 0; i < block_length; i++)
    data[i] = PyInt_AS_LONG(PySequence_Fast_GET_ITEM(block, i));
  
  return data;
}

/**
 * Returns an AmberSystem block as an array of floats.
 *
 * You have to free the thing you get back from this method!
 * Equivalent to: mol.blocks[block_name]
 */
float *get_block_as_float_array(PyObject *mol, const char *block_name, int *size)
{
  PyObject *block;
  if((block = get_block(mol, block_name)) == NULL)
  {
    // fprintf(stderr, "get_block_as_float_array: Couldn't find block %s\n", block_name);
    return NULL;
  }
  
  // Extract the list of doubles from the AmberSystem
  Py_ssize_t block_length = PySequence_Size(block);
  *size = (int) block_length;
  float *data = (float*) malloc(sizeof(float) * (int) block_length);
  for(Py_ssize_t i = 0; i < block_length; i++)
    data[i] = (float) PyFloat_AS_DOUBLE(PySequence_Fast_GET_ITEM(block, i));
  
  return data;
}

/**
 * Equivalent to (mol.x, mol.y, mol.z)
 * You have to free *x, *y, *z!
 */
int get_coordinates(PyObject *mol, float **x, float **y, float **z)
{
  int nx, ny, nz;
  *x = get_float_array(mol, "x", &nx);
  *y = get_float_array(mol, "y", &ny);
  *z = get_float_array(mol, "z", &nz);
  // These coordinate arrays should all be the same length
  if(!(nx == ny && ny == nz))
  {
    fprintf(stderr, "Error: Coordinate arrays are not the same size.\n");
    return -1; 
  }
  return nx;
}

/**
 * Returns the Cartesian distance between two 3D points.
 */
float distance(
  float x0, float y0, float z0,
  float x1, float y1, float z1)
{
  float dx = x1 - x0;
  float dy = y1 - y0;
  float dz = z1 - z0;
  
  return sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Returns the angle (in degrees) formed by these three points,
 * with (x1, y1, z1) forming the vertex of the angle.
 */
float angle(
  float x0, float y0, float z0,
  float x1, float y1, float z1,
  float x2, float y2, float z2)
{
  return 0.0f;
}

/**
 * Calculates the bond energy of the atoms.
 *     E = (1/2) k (l - l0)^2
 * Though for some reason AMBER leaves off the 1/2 coefficient,
 * perhaps because it considers A->B and B->A to be separate, even
 * though they are equal. Civil rights be damned. So, us too.
 */
float calc_bond_energy(
  int *bond_list,
  float *bond_l0,
  float *bond_k,
  float *x,
  float *y,
  float *z,
  int num_atoms,
  int bond_list_size,
  int num_bond_types)
{
  float Ebond = 0.0f;
  
  for(int i = 0; i < bond_list_size; i+=3)
  {
    int ai = bond_list[i] / 3;
    int aj = bond_list[i+1] / 3;
    int bond_id = bond_list[i+2] - 1;

    float l = distance(x[ai], y[ai], z[ai], x[aj], y[aj], z[aj]);
    float dist_diff = l - bond_l0[bond_id];
    float this_Ebond = bond_k[bond_id] * dist_diff * dist_diff;
    Ebond += this_Ebond;
  }
  
  return Ebond;
}

/**
 * Calculates the van der Waals (Lennard-Jones) energy.
 * This is O(N^2) and doesn't bother with cutoffs.
 */
float calc_vdw_energy(
  int *atom_type_index,
  int *nb_parm_index,
  float *lj_a,
  float *lj_b,
  float *x,
  float *y,
  float *z,
  int num_atoms,
  int num_atom_types)
{
  /* From http://ambermd.org/formats.html:
   *  NONBONDED_PARM_INDEX (== ICO) 
   *  provides the index to the nonbon parameter
   *  arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
   *  or 10-12 atoms type interactions are represented.
   *  NOTE: A particular atom type can have either a 10-12
   *  or a 6-12 interaction, but not both.  The index is
   *  calculated as follows:
   *    index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
   *  Where IAC == ATOM_TYPE_INDEX
   *  If index is positive, this is an index into the
   *  6-12 parameter arrays (CN1 and CN2) otherwise it
   *  is an index into the 10-12 parameter arrays (ASOL
   *  and BSOL).
   */
  
  float vdw = 0.0;
  
  for(int ai = 0; ai < num_atoms-1; ai++)
  {
    float xi = x[ai];
    float yi = y[ai];
    float zi = z[ai];

    for(int aj = ai+1; aj < num_atoms; aj++)
    {
      // index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
      int index = nb_parm_index[num_atom_types * (atom_type_index[ai]-1) + atom_type_index[aj]];
      if(index < 0)
      {
        fprintf(stderr, "calc_vdw_energy: Error: Fix me! Negative index %d encountered.\n", index);
        return 0.0f;
      }
      // These A and B values describe the interaction between these
      // two specific atom types.
      float lj_a_ij = lj_a[index];
      float lj_b_ij = lj_b[index];
      float dist = distance(xi, yi, zi, x[aj], y[aj], z[aj]);
      // (A/rij)^12 - 2*(B/rij)^6
      // pow() is really friggin' slow
      //vdw += lj_a_ij / pow(dist, -12) - lj_b_ij / pow(dist, -6);
      float dist3 = dist * dist * dist;
      float dist6 = dist3 * dist3;
      float dist12 = dist6 * dist6;
      vdw += lj_a_ij / dist12 - lj_b_ij / dist6;
    }
  }
  
  return vdw;  
}

/**
 * Calculates the electrostatic energy of the atoms.
 * No fancy stuff like cutoffs and whatnot, just regular old O(N^2).
 */
float calc_elec_energy(
  float *charge,
  float *x,
  float *y,
  float *z,
  int num_atoms)
{
  float eel = 0.0;
  
  for(int ai = 0; ai < num_atoms-1; ai++)
  {
    float xi = x[ai];
    float yi = y[ai];
    float zi = z[ai];
    float qi = charge[ai];
    for(int aj = ai+1; aj < num_atoms; aj++)
    {
      float dist = distance(xi, yi, zi, x[aj], y[aj], z[aj]);
      eel += (qi * charge[aj]) / dist;
    }
  }
  
  return eel;
}
/******** Begin exported functions ********/

static PyObject *mmeval_vdwenergy(PyObject *self, PyObject *args)
{
  // Get the AmberSystem
  PyObject *mol;
  if(!PyArg_ParseTuple(args, "O", &mol))
    return NULL;

  // Extract atom type and L-J parameters from the prmtop
  int header_size = 0;
  int *header = get_block_as_int_array(mol, "POINTERS", &header_size);
  if(header == NULL)
    return NULL;
    
  int nb_parm_index_size = 0;
  int *nb_parm_index = get_block_as_int_array(mol, "NONBONDED_PARM_INDEX", &nb_parm_index_size);
  if(nb_parm_index == NULL)
  {
    free(header);
    return NULL;
  }

  int atom_type_index_size = 0;
  int *atom_type_index = get_block_as_int_array(mol, "ATOM_TYPE_INDEX", &atom_type_index_size);
  if(atom_type_index == NULL)
  {
    free(nb_parm_index);
    free(header);
    return NULL;
  }

  int lj_a_size = 0;
  float *lj_a = get_block_as_float_array(mol, "LENNARD_JONES_ACOEF", &lj_a_size);
  if(lj_a == NULL)
  {
    free(atom_type_index);
    free(nb_parm_index);
    free(header);
    return NULL;
  }

  int lj_b_size = 0;
  float *lj_b = get_block_as_float_array(mol, "LENNARD_JONES_BCOEF", &lj_b_size);
  if(lj_b == NULL)
  {
    free(lj_a);
    free(atom_type_index);
    free(nb_parm_index);
    free(header);
    return NULL;
  }

  // Extract coordinate arrays
  float *x, *y, *z;
  int num_atoms = get_coordinates(mol, &x, &y, &z);
  int num_atom_types = header[1];
  float Evdw = 0.0f;
   
  Evdw = calc_vdw_energy(
    atom_type_index,
    nb_parm_index,
    lj_a,
    lj_b,
    x,
    y,
    z,
    num_atoms,
    num_atom_types);

  free(lj_b);
  free(lj_a);
  free(x);
  free(y);
  free(z);
  free(atom_type_index);
  free(nb_parm_index);
  free(header);
  return Py_BuildValue("f", Evdw);
}

static PyObject *mmeval_elecenergy(PyObject *self, PyObject *args)
{
  // Get the AmberSystem
  PyObject *mol;
  if(!PyArg_ParseTuple(args, "O", &mol))
    return NULL;

  int charge_size = 0;
  float *charge = get_block_as_float_array(mol, "CHARGE", &charge_size);
  if(charge == NULL)
  {
    fprintf(stderr, "Error: Couldn't find CHARGE block.\n");
    return NULL;
  }
  
  // Extract coordinate arrays
  float *x, *y, *z;
  int num_atoms = get_coordinates(mol, &x, &y, &z);
  
  float Eelec = calc_elec_energy(charge, x, y, z, num_atoms);
  free(x);
  free(y);
  free(z);
  return Py_BuildValue("f", Eelec);
}

/**
 * Calculates the bond energy of an AmberSystem.
 */
static PyObject *mmeval_bondenergy(PyObject *self, PyObject *args)
{
  // Get the AmberSystem
  PyObject *mol;
  if(!PyArg_ParseTuple(args, "O", &mol))
    return NULL;
  
  int bond_l0_size = 0, bond_k_size = 0;
  float *bond_l0 = get_block_as_float_array(mol, "BOND_EQUIL_VALUE", &bond_l0_size);
  float *bond_k = get_block_as_float_array(mol, "BOND_FORCE_CONSTANT", &bond_k_size);
  if(bond_l0_size != bond_k_size)
  {
    fprintf(stderr, "Error: Size of BOND_EQUIL_VALUE and BOND_FORCE_CONSTANT blocks is not the same.\n");
    return NULL;
  }
  
  // Extract coordinate arrays
  float *x, *y, *z;
  int num_atoms = get_coordinates(mol, &x, &y, &z);
  
  int bond_list_size = 0;
  int *bond_list = get_block_as_int_array(mol, "BONDS_INC_HYDROGEN", &bond_list_size);
  float Ebond_h = calc_bond_energy(bond_list, bond_l0, bond_k, x, y, z, num_atoms,
    bond_list_size, bond_l0_size);
  free(bond_list);

  bond_list = get_block_as_int_array(mol, "BONDS_WITHOUT_HYDROGEN", &bond_list_size);
  float Ebond_nonh = calc_bond_energy(bond_list, bond_l0, bond_k, x, y, z, num_atoms,
    bond_list_size, bond_l0_size);
  free(bond_list);

  free(x);
  free(y);
  free(z);
  return Py_BuildValue("f", Ebond_h + Ebond_nonh);
}

static PyObject *mmeval_hello(PyObject *self, PyObject *args)
{
  char *foo;

  if (!PyArg_ParseTuple(args, "s", &foo))
    return NULL;
  return Py_BuildValue("i", 42);
}

/* Interface to Python */

static PyMethodDef MMEvalMethods[] =
{
  {"hello", mmeval_hello, METH_VARARGS, "Say hello."},
  {"bondenergy", mmeval_bondenergy, METH_VARARGS, "Evaluate bond energy of an AmberSystem."},
  {"elecenergy", mmeval_elecenergy, METH_VARARGS, "Evaluate electrostatic energy of an AmberSystem."},
  {"vdwenergy", mmeval_vdwenergy, METH_VARARGS, "Evaluate van der Waals energy of an AmberSystem."},
  {NULL, NULL, 0, NULL} // Sentinel
};

PyMODINIT_FUNC initmmeval(void)
{
    (void) Py_InitModule("mmeval", MMEvalMethods);
}
