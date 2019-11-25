//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#include "RotateMesh.h"
#include "core/ClassRegistry.h"

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <vector>
#include <string>
#include <cmath>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, RotateMesh, "rotate_mesh");

RotateMesh::RotateMesh(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    meshPartNames_(),
    meshParts_(),
    rotateVel_(false),
    ndim_(meta_.spatial_dimension())
{
    load(node);
}

void RotateMesh::load(const YAML::Node& node)
{
    const auto& fParts = node["mesh_parts"];
    if (fParts.Type() == YAML::NodeType::Scalar) {
        meshPartNames_.push_back(fParts.as<std::string>());
    } else {
        meshPartNames_ = fParts.as<std::vector<std::string>>();
    }

    angle_ = node["angle"].as<double>();
    angle_ *= std::acos(-1.0) / 180.0;
    axis_ = node["axis"].as<std::vector<double>>();
    origin_ = node["origin"].as<std::vector<double>>();
    if (node["rotate_vel"]) {
        rotateVel_ = node["rotate_vel"].as<bool>();
    }

    ThrowAssert(axis_.size() == 3);
    ThrowAssert(origin_.size() == 3);
}

void RotateMesh::initialize()
{
    for (auto pName: meshPartNames_){
        stk::mesh::Part* part = meta_.get_part(pName);
        if (NULL == part) {
            throw std::runtime_error(
                "RotateMesh: Mesh realm not found in mesh database.");
        } else {
            meshParts_.push_back(part);
        }
    }

    if (rotateVel_) {
        mesh_.add_output_field("velocity");
        mesh_.add_output_field("mesh_displacement");
        mesh_.add_output_field("turbulent_ke");
        mesh_.add_output_field("specific_dissipation_rate");
    }
}

void RotateMesh::run()
{

    bool dowrite = (bulk_.parallel_rank() == 0);
    auto num_steps = mesh_.stkio().get_num_time_steps();
    auto times = mesh_.stkio().get_time_steps();

    auto final_time = times[num_steps - 1];
    std::vector<stk::io::MeshField> missing_fields;
    auto found_time = mesh_.stkio().read_defined_input_fields(final_time, &missing_fields);

    if (missing_fields.size() > 0) {
        if (dowrite) {
            std::cout << "Missing fields in the solution file: " << std::endl;
            for (size_t i=0; i < missing_fields.size(); i++)
                std::cout << "    -" << missing_fields[i].field()->name() << std::endl;
        }
        throw std::runtime_error("ABLStatistics:: missing required fields in database");
    }
    
    if (bulk_.parallel_rank() == 0)
        std::cout << "Rotating mesh " << std::endl;
    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* mesh_disp = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");
    
    stk::mesh::Selector s_part = stk::mesh::selectUnion(meshParts_);
    const stk::mesh::BucketVector& node_buckets = bulk_.get_buckets(
        stk::topology::NODE_RANK, s_part);

    // Calculate the magnitude of the rotation axis vector
    double mag = 0.0;
    for (int i=0; i<3; i++) {
        mag += axis_[i] * axis_[i];
    }
    mag = std::sqrt(mag);

    const int ndim = meta_.spatial_dimension();
    const double cosang = std::cos(0.5*angle_);
    const double sinang = std::sin(0.5*angle_);
    const double q0 = cosang;
    const double q1 = sinang * axis_[0] / mag;
    const double q2 = sinang * axis_[1] / mag;
    const double q3 = sinang * axis_[2] / mag;
    double oldxyz[3] = {0.0, 0.0, 0.0};
    double newxyz[3] = {0.0, 0.0, 0.0};
    double olddisp[3] = {0.0, 0.0, 0.0};
    double newdisp[3] = {0.0, 0.0, 0.0};

    for(auto b: node_buckets) {
        for(size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* xyz = stk::mesh::field_data(*coords, node);
            double* disp = stk::mesh::field_data(*mesh_disp, node);

            for (int i=0; i<ndim; i++) {
                oldxyz[i] = xyz[i];
                olddisp[i] = disp[i];
            }

            const double cx = oldxyz[0] - origin_[0];
            const double cy = oldxyz[1] - origin_[1];
            const double cz = oldxyz[2] - origin_[2];

            const double cdx = olddisp[0] - origin_[0];
            const double cdy = olddisp[1] - origin_[1];
            const double cdz = olddisp[2] - origin_[2];

            newxyz[0] = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * cx +
                2.0 * (q1*q2 - q0*q3) * cy +
                2.0 * (q0*q2 + q1*q3) * cz + origin_[0];

            newxyz[1] = 2.0 * (q1*q2 + q0*q3) * cx +
                (q0*q0 - q1*q1 + q2*q2 - q3*q3) * cy +
                2.0 * (q2*q3 - q0*q1) * cz + origin_[1];

            newxyz[2] = 2.0 * (q1*q3 - q0*q2) * cx +
                2.0 * (q0*q1 + q2*q3) * cy +
                (q0*q0 - q1*q1 - q2*q2 + q3*q3) * cz + origin_[2];

            newdisp[0] = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * cdx +
            2.0 * (q1*q2 - q0*q3) * cdy +
            2.0 * (q0*q2 + q1*q3) * cdz + origin_[0];

            newdisp[1] = 2.0 * (q1*q2 + q0*q3) * cdx +
            (q0*q0 - q1*q1 + q2*q2 - q3*q3) * cdy +
            2.0 * (q2*q3 - q0*q1) * cdz + origin_[1];

            newdisp[2] = 2.0 * (q1*q3 - q0*q2) * cdx +
            2.0 * (q0*q1 + q2*q3) * cdy +
            (q0*q0 - q1*q1 - q2*q2 + q3*q3) * cdz + origin_[2];
            

            for (int i=0; i<ndim; i++) {
                xyz[i] = newxyz[i];
                disp[i] = newdisp[i];                
            }
        }
    }

    if (rotateVel_) {
        VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
            stk::topology::NODE_RANK, "velocity");
        if (velocity != nullptr) {        
            for(auto b: node_buckets) {
                for(size_t in=0; in < b->size(); in++) {
                    auto node = (*b)[in];
                    double* vel = stk::mesh::field_data(*velocity, node);
                    
                    double u = vel[0];
                    double v = vel[1];
                    double w = vel[2];
                    
                    vel[0] = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * u +
                        2.0 * (q1*q2 - q0*q3) * v +
                        2.0 * (q0*q2 + q1*q3) * w;
                    
                    vel[1] = 2.0 * (q1*q2 + q0*q3) * u +
                        (q0*q0 - q1*q1 + q2*q2 - q3*q3) * v +
                        2.0 * (q2*q3 - q0*q1) * w;
                    
                    vel[2] = 2.0 * (q1*q3 - q0*q2) * u +
                        2.0 * (q0*q1 + q2*q3) * v +
                        (q0*q0 - q1*q1 - q2*q2 + q3*q3) * w;
                }
            }
        } else {
            if (bulk_.parallel_rank() == 0) {
                std::cout << "No velocity field found. Skipping " << 
                    "transformation of velocity field" << std::endl;
            }
            rotateVel_ = false;
        }
    }
    
    mesh_.set_write_flag();
}

}  // nalu
}  // sierra
