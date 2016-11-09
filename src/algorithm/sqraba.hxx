//
// Copyright (c) 2016 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_sqraba_hxx__
#define __se3_sqraba_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/sqraba-utils.hxx"

/// @cond DEV

namespace se3
{
  struct SqrAbaForwardStep1 : public fusion::JointVisitor<SqrAbaForwardStep1>
  {
    typedef boost::fusion::vector<const se3::Model &,
    se3::Data &,
    const Eigen::VectorXd &,
    const Eigen::VectorXd &
    > ArgsType;
    
    JOINT_VISITOR_INIT(SqrAbaForwardStep1);
    
    template<typename JointModel>
    static void algo(const se3::JointModelBase<JointModel> & jmodel,
                     se3::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const se3::Model & model,
                     se3::Data & data,
                     const Eigen::VectorXd & q,
                     const Eigen::VectorXd & v)
    {
      const Model::JointIndex & i = jmodel.id();
      jmodel.calc(jdata.derived(),q,v);
      
      const Model::Index & parent = model.parents[i];
      data.liMi[i] = model.jointPlacements[i] * jdata.M();
      
      data.v[i] = jdata.v();
      
      if (parent>0)
      {
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
        data.v[i] += data.liMi[i].actInv(data.v[parent]);
      }
      else
        data.oMi[i] = data.liMi[i];
      
      data.a[i] = jdata.c() + (data.v[i] ^ jdata.v());
      
      data.Yaba[i] = model.inertias[i].matrix();
      data.f[i] = model.inertias[i].vxiv(data.v[i]); // -f_ext
    }
    
  };
  
  struct SqrAbaBackwardStep : public fusion::JointVisitor<SqrAbaBackwardStep>
  {
    typedef boost::fusion::vector<const Model &,
    Data &> ArgsType;
    
    JOINT_VISITOR_INIT(SqrAbaBackwardStep);
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data)
    {
      const Model::JointIndex & i = jmodel.id();
      const Model::Index & parent  = model.parents[i];
      Inertia::Matrix6 & Ia = data.Yaba[i];
    }
    
  };
  
  struct SqrAbaForwardStep2 : public fusion::JointVisitor<SqrAbaForwardStep2>
  {
    typedef boost::fusion::vector<const se3::Model &,
    se3::Data &
    > ArgsType;
    
    JOINT_VISITOR_INIT(SqrAbaForwardStep2);
    
    template<typename JointModel>
    static void algo(const se3::JointModelBase<JointModel> & jmodel,
                     se3::JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const se3::Model & model,
                     se3::Data & data)
    {
      const Model::JointIndex & i = jmodel.id();
      const Model::Index & parent = model.parents[i];
    }
    
  };
  
  inline const Eigen::VectorXd &
  sqraba(const Model & model,
         Data & data,
         const Eigen::VectorXd & q,
         const Eigen::VectorXd & v,
         const Eigen::VectorXd & tau)
  {
    data.v[0].setZero();
    data.a[0] = -model.gravity;
    data.u = tau;
    
    for(Model::Index i=1;i<(Model::Index)model.njoints;++i)
    {
      SqrAbaForwardStep1::run(model.joints[i],data.joints[i],
                              SqrAbaForwardStep1::ArgsType(model,data,q,v));
    }
    
    for( Model::Index i=(Model::Index)model.njoints-1;i>0;--i )
    {
      //SqrAbaBackwardStep::run(model.joints[i],data.joints[i],
      //SqrAbaBackwardStep::ArgsType(model,data));
    }
    
    for(Model::Index i=1;i<(Model::Index)model.njoints;++i)
    {
      //SqrAbaForwardStep2::run(model.joints[i],data.joints[i],
      //                      SqrAbaForwardStep2::ArgsType(model,data));
    }
    
    return data.ddq;
  }

} // namespace se3

/// @endcond

#endif // ifndef __se3_sqraba_hxx__
