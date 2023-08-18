//
// Created by joaovictor on 17/08/23.
//

#ifndef GERVLIB_PAGEMANAGER_H
#define GERVLIB_PAGEMANAGER_H

#include "Page.h"
#include "Utils.h"
#include "Configure.h"
#include <map>
#include <vector>
#include <fstream>

namespace gervLib::memory
{

    template <typename O>
    class PageManager : public serialize::Serialize{

    private:
        O pageSize;
        std::map<O, std::vector<std::unique_ptr<Page<O>>>> pages; // node id -> [(page id, start, end), ...]
        std::string pagePrefix, folder;
        std::unique_ptr<Page<O>> currentPage;

    public:
        explicit PageManager(std::string _pagePrefix, std::string _folder, O _pageSize = 0) : pageSize(_pageSize), pagePrefix(std::move(_pagePrefix)), folder(std::move(_folder))
        {

            utils::createFolderIfNotExists(folder);

            if (pageSize != 0)
            {
                currentPage = std::make_unique<Page<O>>();
                currentPage->setPath(gervLib::utils::generateFileByPrefix(folder, pagePrefix, true, ".bin"));
            }
            else
            {
                currentPage = nullptr;
            }

        }

        PageManager()
        {
            currentPage = nullptr;
        }

        virtual ~PageManager()
        {

            for (auto &page : pages)
            {
                for (auto &p : page.second)
                {
                    p.reset();
                }
            }

            pages.clear();

        }

        //getSerializedSize
        //serialize
        //deserialize

        std::unique_ptr<u_char[]> serialize() override
        {

            std::unique_ptr<u_char[]> data(new u_char[getSerializedSize()]);
            size_t offset = 0, size, pages_len = pages.size();

            std::memcpy(data.get() + offset, &pageSize, sizeof(O));
            offset += sizeof(O);

            if (currentPage == nullptr)
            {
                size = 0;
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
            }
            else
            {
                size = currentPage->getSerializedSize();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);
                std::unique_ptr<u_char[]> data2 = currentPage->serialize();
                std::memcpy(data.get() + offset, data2.get(), size);
                offset += size;
                data2.reset();
            }

            size = pagePrefix.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, pagePrefix.c_str(), size);
            offset += size;

            size = folder.size();
            std::memcpy(data.get() + offset, &size, sizeof(size_t));
            offset += sizeof(size_t);

            std::memcpy(data.get() + offset, folder.c_str(), size);
            offset += size;

            std::memcpy(data.get() + offset, &pages_len, sizeof(size_t));
            offset += sizeof(size_t);

            for (auto &page : pages)
            {

                std::memcpy(data.get() + offset, &page.first, sizeof(O));
                offset += sizeof(O);

                size = page.second.size();
                std::memcpy(data.get() + offset, &size, sizeof(size_t));
                offset += sizeof(size_t);

                for (auto &p : page.second)
                {

                    size = p->getSerializedSize();
                    std::memcpy(data.get() + offset, &size, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> data2 = p->serialize();
                    std::memcpy(data.get() + offset, data2.get(), size);
                    offset += size;
                    data2.reset();

                }

            }

            return data;

        }

        void deserialize(std::unique_ptr<u_char[]> _data) override
        {

            size_t offset = 0, size, pages_len;

            std::memcpy(&pageSize, _data.get() + offset, sizeof(O));
            offset += sizeof(O);

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            if (size == 0)
            {
                currentPage = nullptr;
            }
            else
            {
                currentPage = std::make_unique<Page<O>>();
                std::unique_ptr<u_char[]> data2(new u_char[size]);
                std::memcpy(data2.get(), _data.get() + offset, size);
                offset += size;
                currentPage->deserialize(std::move(data2));
            }

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            pagePrefix.resize(size);
            std::memcpy(&pagePrefix[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&size, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            folder.resize(size);
            std::memcpy(&folder[0], _data.get() + offset, size);
            offset += size;

            std::memcpy(&pages_len, _data.get() + offset, sizeof(size_t));
            offset += sizeof(size_t);

            for (size_t i = 0; i < pages_len; i++)
            {

                O nodeID;
                std::memcpy(&nodeID, _data.get() + offset, sizeof(O));
                offset += sizeof(O);

                std::vector<std::unique_ptr<Page<O>>> v;
                size_t v_len;
                std::memcpy(&v_len, _data.get() + offset, sizeof(size_t));
                offset += sizeof(size_t);

                for (size_t j = 0; j < v_len; j++)
                {

                    size_t size2;
                    std::memcpy(&size2, _data.get() + offset, sizeof(size_t));
                    offset += sizeof(size_t);

                    std::unique_ptr<u_char[]> data2(new u_char[size2]);
                    std::memcpy(data2.get(), _data.get() + offset, size2);
                    offset += size2;

                    std::unique_ptr<Page<O>> p = std::make_unique<Page<O>>();
                    p->deserialize(std::move(data2));
                    v.push_back(std::move(p));

                }

                pages[nodeID] = std::move(v);

            }

            _data.reset();

        }

        size_t getSerializedSize() override
        {

            size_t ans = sizeof(size_t);

            if (currentPage != nullptr)
                ans += sizeof(size_t) + currentPage->getSerializedSize();
            else
                ans += sizeof(size_t);

            ans += sizeof(size_t) + pagePrefix.size();
            ans += sizeof(size_t) + folder.size();
            ans += sizeof(size_t);

            for (auto &page : pages)
            {

                ans += sizeof(O);
                ans += sizeof(size_t);

                for (auto &p : page.second)
                {

                    ans += sizeof(size_t);
                    ans += p->getSerializedSize();

                }

            }

            return ans;

        }

        bool isEqual(const PageManager &other) const
        {

            if (pageSize != other.pageSize)
                return false;

            if (pagePrefix != other.pagePrefix)
                return false;

            if (folder != other.folder)
                return false;

            if (pages.size() != other.pages.size())
                return false;

            for (auto &page : pages)
            {

                if (other.pages.find(page.first) == other.pages.end())
                    return false;

                if (page.second.size() != other.pages.at(page.first).size())
                    return false;

                for (size_t i = 0; i < page.second.size(); i++)
                {

                    if (*page.second[i] != *other.pages.at(page.first)[i])
                        return false;

                }

            }

            return true;

        }

        friend std::ostream &operator<<(std::ostream &os, const PageManager &manager)
        {

            os << "PageManager: " << std::endl;
            os << "Page size: " << manager.pageSize << std::endl;
            os << "Page prefix: " << manager.pagePrefix << std::endl;
            os << "Folder: " << manager.folder << std::endl;

            if (manager.currentPage == nullptr)
                os << "Current page: NULL" << std::endl;
            else
                os << "Current page: " << *manager.currentPage << std::endl;

            os << "Pages: " << std::endl;

            for (auto &page : manager.pages)
            {
                os << "Node " << page.first << ": " << std::endl;

                for (auto &p : page.second)
                {
                    os << *p << std::endl;
                }

            }

            return os;

        }

        O getPageSize() const
        {
            return pageSize;
        }

        std::vector<std::unique_ptr<Page<O>>> &getPages(O nodeId)
        {
            return pages[nodeId];
        }

        void save(O nodeID, std::unique_ptr<u_char[]> data, O size)
        {

            if (pages.find(nodeID) != pages.end())
                return;

            if (pageSize == 0)
            {

                std::string path = utils::generateFileByPrefix(folder, pagePrefix, true, ".bin");
                std::ofstream file(path, std::ios::binary);
                configure::IOWrite++;

                if (file.is_open())
                {
                    file.write((char *) data.get(), size);
                    file.close();
                    pages[nodeID].push_back(std::make_unique<Page<O>>(0, size, path));
                }
                else
                    throw std::runtime_error("PageManager::save(): Could not open file " + path);

            }
            else
            {

                size_t start, end, currentCopy = 0, cpy;
                std::string path;

                while (true)
                {

                    start = (currentPage->empty()) ? 0 : currentPage->getEnd() + 1;
                    end = start + (size - currentCopy) - 1;
                    cpy = (end > pageSize) ? pageSize - start : size - currentCopy;
                    end = std::min(end, pageSize - 1);
                    path = currentPage->getPath();

                    std::ofstream file(path, std::ios::binary | std::ios::app);
                    configure::IOWrite++;

                    if (file.is_open())
                    {
                        file.seekp(start);
                        file.write((char *) data.get() + currentCopy, cpy);
                        file.close();
                        pages[nodeID].push_back(std::make_unique<Page<O>>(start, end, path));
                        currentCopy += cpy;
                        currentPage->setEnd(end);
                    }
                    else
                        throw std::runtime_error("PageManager::save(): Could not open file " + path);

                    if (currentPage->getEnd() == (pageSize - 1))
                    {
                        currentPage = std::make_unique<Page<O>>();
                        currentPage->setPath(utils::generateFileByPrefix(folder, pagePrefix, true, ".bin"));
                    }

                    if (currentCopy == size)
                        break;

                }

            }

            data.reset();

        }

        std::unique_ptr<u_char[]> load(O nodeID)
        {

            O total = 0;

            if (pages.find(nodeID) == pages.end())
                throw std::runtime_error("PageManager::load(): Node " + std::to_string(nodeID) + " not found");

            for (auto &page : pages[nodeID])
            {
                total += page->getEnd() - page->getStart() + 1;
            }

            std::unique_ptr<u_char[]> data(new u_char[total]);
            O offset = 0, size;

            for (auto &page : pages[nodeID])
            {

                std::ifstream file(page->getPath(), std::ios::binary);
                configure::IORead++;

                if (file.is_open())
                {
                    file.seekg(page->getStart());
                    size = page->getEnd() - page->getStart() + 1;
                    file.read((char *) data.get() + offset, size);
                    file.close();
                    offset += size;
                }
                else
                    throw std::runtime_error("PageManager::load(): Could not open file " + page->getPath());

            }

            return data;

        }




    };

}





#endif //GERVLIB_PAGEMANAGER_H
