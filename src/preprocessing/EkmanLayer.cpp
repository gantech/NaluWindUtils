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


#include "EkmanLayer.h"
#include "core/LinearInterpolation.h"

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, EkmanLayer, "init_ekman_layer");

EkmanLayer::EkmanLayer(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh),
    meta_(mesh.meta()),
    bulk_(mesh.bulk()),
    fluid_parts_(0),
    ndim_(meta_.spatial_dimension()),
    doVelocity_(false),
    doVelocityPerturb_(false)
{
    load(node);
}

void EkmanLayer::load(const YAML::Node& ekmanLayer)
{
    auto fluid_partnames = ekmanLayer["fluid_parts"].as<std::vector<std::string>>();

    if (ekmanLayer["velocity"]) {
        doVelocity_ = true;
        load_velocity_info(ekmanLayer["velocity"]);
    }

    fluid_parts_.resize(fluid_partnames.size());

    for(size_t i=0; i<fluid_partnames.size(); i++) {
        stk::mesh::Part* part = meta_.get_part(fluid_partnames[i]);
        if (NULL == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    }
}

void EkmanLayer::initialize()
{
    if (doVelocity_) {
        VectorFieldType& velocity = meta_.declare_field<VectorFieldType>(
            stk::topology::NODE_RANK, "velocity");
        for(auto part: fluid_parts_) {
            stk::mesh::put_field(velocity, *part);
        }
        mesh_.add_output_field("velocity");
    }
}

void EkmanLayer::run()
{
    if (bulk_.parallel_rank() == 0)
        std::cout << "Generating Ekman Layer fields" << std::endl;
    if (doVelocity_) init_velocity_field();

}

void EkmanLayer::load_velocity_info(const YAML::Node& ekmanLayer)
{

    G_ = ekmanLayer["G"].as<double>();

    alpha_ = ekmanLayer["alpha"].as<double>();

    D_ = ekmanLayer["D"].as<double>();

    ky_ = ekmanLayer["ky"].as<double>();

    if (ekmanLayer["perturb"]) {

        doVelocityPerturb_ = true ;
        
        avxHeights_ = ekmanLayer["avxHeights"].as<std::vector<double>>();
        auto navxHeights = avxHeights_.size();
        std::cout << "navxHeights = " << navxHeights  << std::endl << std::flush;

        auto avxData = ekmanLayer["avx"].as<std::vector<std::vector<double>>>();
        ThrowAssertMsg(
            (navxHeights == avxData.size()),
            "EkmanLayer: Mismatch between sizes of heights and velocity amplitude in X provided "
            "for initializing EKMANLAYER fields. Check input file.");

        ThrowAssertMsg(
            (2 == avxData.at(0).size()),
            "EkmanLayer : Velocity X perturbation amplitude should be complex number with 2 components");

        avx_.resize(navxHeights);
        for (int i=0; i < navxHeights; i++) {
            avx_[i] = avxData[i][0] + 1i * avxData[i][1];
        }

        avyHeights_ = ekmanLayer["avyHeights"].as<std::vector<double>>();
        auto navyHeights = avyHeights_.size();

        auto avyData = ekmanLayer["avy"].as<std::vector<std::vector<double>>>();
        ThrowAssertMsg(
            (navyHeights == avyData.size()),
            "EkmanLayer: Mismatch between sizes of heights and velocity amplitude in Y provided "
            "for initializing EKMANLAYER fields. Check input file.");

        ThrowAssertMsg(
            (2 == avyData.at(0).size()),
            "EkmanLayer : Velocity Y perturbation amplitude should be complex number with 2 components");

        avy_.resize(navyHeights);
        for (int i=0; i < navyHeights; i++) {
            avy_[i] = avyData[i][0] + 1i * avyData[i][1];
        }
        
        avzHeights_ = ekmanLayer["avzHeights"].as<std::vector<double>>();
        auto navzHeights = avzHeights_.size();

        auto avzData = ekmanLayer["avz"].as<std::vector<std::vector<double>>>();
        ThrowAssertMsg(
            (navzHeights == avzData.size()),
            "EkmanLayer: Mismatch between sizes of heights and velocity amplitude in Z provided "
            "for initializing EKMANLAYER fields. Check input file.");
        
        ThrowAssertMsg(
            (2 == avzData.at(0).size()),
            "EkmanLayer : Velocity Z perturbation amplitude should be complex number with 2 components");

        avz_.resize(navzHeights);
        for (int i=0; i < navzHeights; i++) {
            avz_[i] = avzData[i][0] + 1i * avzData[i][1];
        }
    }
        
}

void EkmanLayer::init_velocity_field()
{
    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    const stk::mesh::BucketVector& fluid_bkts = bulk_.get_buckets(
        stk::topology::NODE_RANK, fluid_union);

    VectorFieldType* coords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* velocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "velocity");

    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* vel = stk::mesh::field_data(*velocity, fbkt);

        for (size_t in=0; in < fbkt.size(); in++) {
            const double xh = xyz[in*ndim_ + 0];
            const double yh = xyz[in*ndim_ + 1];
            const double zh = xyz[in*ndim_ + 2];

            vel[in * ndim_ + 0] = G_ * ( cos(alpha_) - exp(-zh/D_) * cos(zh/D_ - alpha_) ) ;
            vel[in * ndim_ + 1] = G_ * ( sin(alpha_) + exp(-zh/D_) * sin(zh/D_ - alpha_) ) ;
            vel[in * ndim_ + 2] = 0.0 ;

            if (doVelocityPerturb_) {

                std::complex<double> avx, avy, avz;
                
                utils::linear_interp(
                    avxHeights_, avx_, zh, avx);
                utils::linear_interp(
                    avyHeights_, avy_, zh, avy);
                utils::linear_interp(
                    avzHeights_, avz_, zh, avz);

                std::complex<double> expiky = cos(ky_ * yh) + 1i * sin(ky_ * yh);

                vel[in * ndim_ + 0] += ( avx * expiky ).real();
                vel[in * ndim_ + 1] += ( avy * expiky ).real();
                vel[in * ndim_ + 2] += ( avz * expiky ).real();
                
            }
        }
    }
}

} // nalu
} // sierra
