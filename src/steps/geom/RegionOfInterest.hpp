#pragma once

#include <map>

#include "util/vocabulary.hpp"

namespace steps::tetmesh {

enum ROIType { ROI_VERTEX, ROI_TRI, ROI_TET, ROI_UNDEFINED = 99 };

template <typename DATA_TYPE>
struct ROIBaseTypeTraits {
    using element_type = DATA_TYPE;
    using data_type = std::vector<element_type>;
    using roi_map_type = std::map<std::string, std::vector<element_type>>;
    using insert_return_type = std::pair<typename roi_map_type::iterator, bool>;
};

template <ROIType>
struct ROITypeTraits {};

template <>
struct ROITypeTraits<ROI_VERTEX>: public ROIBaseTypeTraits<vertex_id_t> {
    template <typename ROI>
    static const roi_map_type& get_container(const ROI& roi) {
        return roi.vertices_roi;
    }
    template <typename ROI>
    static roi_map_type& get_container(ROI& roi) {
        return roi.vertices_roi;
    }
};

template <>
struct ROITypeTraits<ROI_TRI>: public ROIBaseTypeTraits<triangle_global_id> {
    template <typename ROI>
    static const roi_map_type& get_container(const ROI& roi) {
        return roi.tris_roi;
    }
    template <typename ROI>
    static roi_map_type& get_container(ROI& roi) {
        return roi.tris_roi;
    }
};

template <>
struct ROITypeTraits<ROI_TET>: public ROIBaseTypeTraits<tetrahedron_global_id> {
    template <typename ROI>
    static const roi_map_type& get_container(const ROI& roi) {
        return roi.tets_roi;
    }
    template <typename ROI>
    static roi_map_type& get_container(ROI& roi) {
        return roi.tets_roi;
    }
};

struct RegionOfInterest {
    template <ROIType value>
    typename ROITypeTraits<value>::roi_map_type::const_iterator get(const std::string& id,
                                                                    unsigned int count = 0,
                                                                    bool warning = true) const;

    template <ROIType value>
    inline typename ROITypeTraits<value>::insert_return_type insert(
        const std::string& id,
        const typename ROITypeTraits<value>::data_type& data) {
        return ROITypeTraits<value>::get_container(*this).emplace(id, data);
    }

    template <ROIType value>
    inline bool replace(const std::string& id,
                        const std::vector<typename ROITypeTraits<value>::element_type>& data) {
        auto it = ROITypeTraits<value>::get_container(*this).find(id);
        if (it == end<value>()) {
            return false;
        }
        it->second.assign(data.begin(), data.end());
        return true;
    }

    template <ROIType value>
    inline std::size_t erase(const std::string& id) {
        return ROITypeTraits<value>::get_container(*this).erase(id);
    }

    template <ROIType value>
    inline typename ROITypeTraits<value>::roi_map_type::const_iterator end() const {
        return ROITypeTraits<value>::get_container(*this).end();
    }

    std::size_t size() const {
        return tets_roi.size() + tris_roi.size() + vertices_roi.size();
    }

    void ids(std::vector<std::string>& ids) const {
        ids.reserve(size());
        for (auto const& it: tets_roi) {
            ids.push_back(it.first);
        }
        for (auto const& it: tris_roi) {
            ids.push_back(it.first);
        }
        for (auto const& it: vertices_roi) {
            ids.push_back(it.first);
        }
    }

    ROITypeTraits<ROI_TET>::roi_map_type tets_roi;
    ROITypeTraits<ROI_TRI>::roi_map_type tris_roi;
    ROITypeTraits<ROI_VERTEX>::roi_map_type vertices_roi;
};

// explicit template instantiation declarations
extern template ROITypeTraits<ROI_VERTEX>::roi_map_type::const_iterator
RegionOfInterest::get<ROI_VERTEX>(const std::string& id, unsigned int count, bool warning) const;
extern template ROITypeTraits<ROI_TRI>::roi_map_type::const_iterator
RegionOfInterest::get<ROI_TRI>(const std::string& id, unsigned int count, bool warning) const;
extern template ROITypeTraits<ROI_TET>::roi_map_type::const_iterator
RegionOfInterest::get<ROI_TET>(const std::string& id, unsigned int count, bool warning) const;

}  // namespace steps::tetmesh
