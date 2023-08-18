//
// Created by joaovictor on 17/08/23.
//

#ifndef GERVLIB_PAGE_H
#define GERVLIB_PAGE_H

#include <cstring>
#include "Serialize.h"


namespace gervLib::memory {

    template<typename O>
    class Page : public gervLib::serialize::Serialize {
    private:
        O start, end; //last valid character
        std::string path;

    public:
        Page() : start(0), end(0) {}

        Page(O _start, O _end, std::string _path) : start(_start), end(_end), path(std::move(_path)) {}

        Page(const Page &other) : start(other.start), end(other.end), path(other.path) {}

        virtual ~Page() = default;

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size;

            std::memcpy(data.get() + offset, &start, sizeof(O));
            offset += sizeof(O);

            std::memcpy(data.get() + offset, &end, sizeof(O));
            offset += sizeof(O);

            size = path.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, path.c_str(), size);
            offset += size;

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, size;

            std::memcpy(&start, _data.get() + offset, sizeof(O));
            offset += sizeof(O);

            std::memcpy(&end, _data.get() + offset, sizeof(O));
            offset += sizeof(O);

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            path.resize(size);
            std::memcpy(&path[0], _data.get() + offset, size);
            offset += size;

            _data.reset();

        }

        size_t getSerializedSize() override
        {
            return sizeof(O) * 3 + path.size();
        }

        bool empty() {

            return end == start;

        }

        bool isEqual(Page<O> &other) {

            return ((start == other.start) && (end == other.end) && (path == other.path));

        }

        bool operator==(Page<O> &other) {

            return isEqual(other);

        }

        bool operator!=(Page<O> &other) {

            return !isEqual(other);

        }

        O getStart() const {

            return start;

        }

        O getEnd() const {

            return end;

        }

        [[nodiscard]] std::string getPath() const {

            return path;

        }

        void getStart(O _start) {

            start = _start;

        }

        void setEnd(O _end) {

            end = _end;

        }

        void setPath(std::string _path) {

            path = std::move(_path);

        }

        friend std::ostream &operator<<(std::ostream &os, const Page &page) {

            os << "Page_t: start = " << page.start << ", end = " << page.end << ", path = " << page.path;
            return os;

        }

    };

}

#endif //GERVLIB_PAGE_H
