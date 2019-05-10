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

#ifndef ROTATEMESH_H
#define ROTATEMESH_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/** Rotate vector fields on a mesh. Optionally add a `delta` before rotation.
 *
 *  ```
 *  rotate_fields:
 *    mesh_parts:
 *      - unspecified-2-hex
 *
 *    field_names:
 *      - name: velocity
 *        delta: [-8.0, 0.0, 0.0]
 *
 *    angle: 45.0
 *    axis: [0.0, 0.0, 1.0]
 *  ```
 */
class RotateFields: public PreProcessingTask
{
public:
    RotateFields(CFDMesh&, const YAML::Node&);

    virtual ~RotateFields() {}

    virtual void initialize();

    virtual void run();

private:
    RotateFields() = delete;
    RotateFields(const RotateFields&) = delete;

    void load(const YAML::Node&);

    void populate_solution();

    //! STK Metadata object
    stk::mesh::MetaData& meta_;

    //! STK Bulkdata object
    stk::mesh::BulkData& bulk_;

    //! Part names of the mesh that needs to be rotated
    std::vector<std::string> meshPartNames_;

    //! Parts of the mesh that need to be rotated
    stk::mesh::PartVector meshParts_;

    //! Names of the fields that need to be rotated
    std::vector<std::string> fieldNames_;

    //! Delta to be added for each field before rotation
    std::vector<std::array<double,3>> fieldDelta_;
    
    //! Angle of rotation
    double angle_;

    //! Axis around which the rotation is performed
    std::vector<double> axis_;

    //! Dimensionality of the mesh
    int ndim_;
};

} // namespace nalu
} // namespace sierra

#endif /* ROTATEMESH_H */
