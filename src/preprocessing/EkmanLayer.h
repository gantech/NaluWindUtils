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

#ifndef EKMANLAYER_H
#define EKMANLAYER_H

#include "PreProcessingTask.h"
#include <complex>

namespace sierra {
namespace nalu {

/**
 * Initialize velocity field to the Ekman layer profile with some perturbation waves in the y direction
 *
 * This task is activated by using the `init_ekman_layer` task in the
 * preprocessing input file. It requires a section `init_ekman_layer` in the
 * `nalu_preprocess` section with the following parameters:
 * 
 * ```
 * U(y,z) = G_ * ( cos(alpha_) - exp(-z/D_) * cos(z/D_ - alpha_) ) + real( avx_(z) * exp(i * k * y ) )
 * V(y,z) = G_ * ( sin(alpha_) + exp(-z/D_) * sin(z/D_ - alpha_) ) + real( avy_(z) * exp(i * k * y ) )
 * W(y,z) = 0.0 + real( avz_(z) * exp(i * k * y) )
 * ```
 * 
 * + G_ - Geostrophic wind (m/s)
 * + D_ - Ekman layer height (m)
 * + alpha_ - Angle of Geostrophic wind with respect to the 'x' direction (Ideally, it should be the east direction. In the future, I may read in east and north direction vectors like the rest of Nalu) (radians)
 * + ky - Wavenumber of inital wave in the 'y' direction (1/m)
 * + avx_, avy_ and avz_ - Complex amplitudes of initial fluctuations as a function of height. (m/s)
 * + avxHeights_, avyHeights_, avzHeights_ - The heights at which the initial fluctuation amplitudes are prescribed. (m)
 *
 *  ```
 *  init_ekman_layer:
 *  fluid_parts: [Unspecified-2-HEX]
 * 
 *  velocity:
 *    G: 1.0
 *    alpha: 0.0
 *    D: 1.0
 *    ky: 0.5
 * 
 *    perturb: True
 * 
 *      avxHeights: [0.0, 10.0, 30.0, 70.0, 100.0, 650.0, 10000.0]
 *      avx:
 *        - [ 0.0, 0.0]
 *        - [4.81947, -4.81947]
 *        - [5.63845, -5.63845]
 *        - [6.36396, -6.36396]
 *        - [6.69663, -6.69663]
 *        - [8.74957, -8.74957]
 *        - [8.74957, -8.74957]
 * 
 *      avyHeights: [0.0, 10.0, 30.0, 70.0, 100.0, 650.0, 10000.0]
 *      avy:
 *        - [ 0.0, 0.0]
 *        - [4.81947, -4.81947]
 *        - [5.63845, -5.63845]
 *        - [6.36396, -6.36396]
 *        - [6.69663, -6.69663]
 *        - [8.74957, -8.74957]
 *        - [8.74957, -8.74957]
 * 
 *      avzHeights: [0.0, 10.0, 30.0, 70.0, 100.0, 650.0, 10000.0]
 *      avz:
 *        - [ 0.0, 0.0]
 *        - [4.81947, -4.81947]
 *        - [5.63845, -5.63845]
 *        - [6.36396, -6.36396]
 *        - [6.69663, -6.69663]
 *        - [8.74957, -8.74957]
 *        - [8.74957, -8.74957]
 *  ```
 *
 * The sections `perturb` is optional, allowing the user to
 * initialize only the laminar velocity profile or including perturbations as desired. 
 * The heights are in meters, the complex amplitude of the velocity is in m/s.
 */
    
class EkmanLayer: public PreProcessingTask
{
public:
    template<typename T>
    using ArrayComplex = std::vector<std::complex<T>>;

    EkmanLayer(CFDMesh&, const YAML::Node&);

    virtual ~EkmanLayer() {}

    //! Declare velocity and temperature fields and register them for output
    void initialize();

    //! Initialize the velocity and/or temperature fields by linear interpolation
    void run();

private:
    EkmanLayer() = delete;
    EkmanLayer(const EkmanLayer&) = delete;

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    //! Helper function to parse and initialize velocity inputs
    void load_velocity_info(const YAML::Node&);

    //! Initialize the velocity field through linear interpolation
    void init_velocity_field();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Parts of the fluid mesh where velocity/temperature is initialized
    stk::mesh::PartVector fluid_parts_;

    //! Geostropic wind
    double G_;

    //! Angle of geostrophic wind w.r.t east direction
    double alpha_;

    //! Ekman Layer height
    double D_;

    //! Wavenumber of modes to create along the y direction
    double ky_;
    
    //! List of heights where complex amplitude of velocity X perturbations is defined
    std::vector<double> avxHeights_;

    //! List of velocity X perturbation amplitudes (complex) at the user-defined heights
    ArrayComplex<double> avx_;

    //! List of heights where complex amplitude of velocity Y perturbations is defined
    std::vector<double> avyHeights_;

    //! List of velocity Y perturbation amplitudes (complex) at the user-defined heights
    ArrayComplex<double> avy_;

    //! List of heights where complex amplitude of velocity Z perturbations is defined
    std::vector<double> avzHeights_;

    //! List of velocity X perturbation amplitudes (complex) at the user-defined heights
    ArrayComplex<double> avz_;

    //! Dimensionality of the mesh
    int ndim_;

    //! Flag indicating whether velocity is initialized
    bool doVelocity_;

    //! Flag indicating whether perturbations in velocity is initialized
    bool doVelocityPerturb_;

};

}
}

#endif /* EKMANLAYER_H */
