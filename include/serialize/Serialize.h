//
// Created by joaovictor on 15/08/23.
//

#ifndef GERVLIB_SERIALIZE_H
#define GERVLIB_SERIALIZE_H

#include <string>
#include <memory>

namespace gervLib::serialize
{

    class Serialize
    {

    public:
        virtual std::unique_ptr<u_char[]> serialize() = 0;
        virtual void deserialize(std::unique_ptr<u_char[]> _data) = 0;
        virtual size_t getSerializedSize() = 0;

    };

}

#endif //GERVLIB_SERIALIZE_H
