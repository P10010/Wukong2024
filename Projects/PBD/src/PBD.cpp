#include "../include/PBD.h"
#include <igl/readOBJ.h>

#include "../include/Util.h"
void PBD::initializeFromFile(const std::string& filename) {
  MatrixXT V;
  MatrixXi F;
  igl::readOBJ(filename, V, F);

  //these functions copy the V, F, that is, what was read on the obj, to the internal representation atRest and faces
  //they were copied from DiscreteShell and also "fatten" it, but that is probably unnecessary for us.
  //you could set the second parameter to 1 rather than 3 if you didn't want to fatten (but at that point it might be worth it creating another function in Util.h)
  //you would probably need to change app.initializeScene too
  //maybe fatten is actually necessary to give some "thickness" to the surface mesh? idk. Changing it breaks something else, I'll leave it as is for now
  iglMatrixFatten<T, 3>(V, atRest);
  iglMatrixFatten<int, 3>(F, faces);
  currentV=atRest;
  int nRows=currentV.size();
  w.resize(nRows);
  v.resize(nRows);
  //todo put 1/mass, at the moment I set the mass to be 1 for every vertex
  for(int i=0;i<nRows;i++)
    w(i)=1;

  for(int i=0;i<nRows;i++)
    v(i)=0;
}

bool PBD::advanceOneStep(int step)
{
  int nRows=currentV.size();
  for(int i=0;i<nRows;i++)
    v(i)+=dt*g;

  //todo line 6, damp velocities
  for(int i=0;i<nRows;i++)
    //todo this should update p(i), not currentV(i)
    currentV(i)+=dt*v(i);

  if(step==10)
    return true;
  return false;

}
