/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "GeneralVectorPostprocessor.h"
#include "MyTRIMRasterizer.h"

class PKAList;

template <>
InputParameters validParams<PKAList>();

/**
 * Outputs the list of MyTRIM defects comopiled with the the MyTRIMDiracRunner
 */
class PKAList : public GeneralVectorPostprocessor
{
public:
  PKAList(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();

protected:
  const MyTRIMRasterizer & _rasterizer;

  /// primary knock-on atom (PKA) list
  const std::vector<MyTRIM_NS::IonBase> & _pka_list;

  ///@{ PKA data to output
  VectorPostprocessorValue & _x;
  VectorPostprocessorValue & _y;
  VectorPostprocessorValue & _z;
  VectorPostprocessorValue & _seed;
  VectorPostprocessorValue & _m;
  VectorPostprocessorValue & _Z;
  ///@}
};

