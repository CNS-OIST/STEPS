#pragma once

#include <queue>
#include <unordered_map>
#include <vector>

namespace steps::util {

/**
 * \brief Priority queue
 *
 * Convenience class. Conceived to do the bookkeeping during some tetrahedron searches.
 * It is generic enough to be used for any A* search. The main capabilities of this
 * class are:
 *
 * - mark down what was investigated
 * - add keys as you please: they will be added only once and investigated only once.
 * - next: provides the next element in the queue based on the distance metric added with
 *   push. The value with the smallest d comes first.
 * - ability to mark any T as investigated and skip it if necessary if it appears in the queue.
 * - keeps track of the next value of T that was not yet added to the queue.
 *
 * The standard way to use this class is this:
 * PQueue pq;
 * pq.push(d, k);
 * while(!pq.empty()) {
 *    auto n = pq.next();
 *    // operate on n, add to the queue, and mark stuff as investigated as you please
 * }
 *
 * There is no size() because it may be unprecise.
 */
template <class D, class T>
class PQueue {
  public:
    /// @brief Status of the particular key. The key is often a tet_id.
    enum class StatusType { UNQUEUED, QUEUED, PROCESSED };

    /// @brief We need to provide the min unprocessed value if different from 0.
    PQueue(const T& min_unqueued_value = T(0))
        : min_unqueued_value_(min_unqueued_value) {}

    /// @brief Retrieve status.
    StatusType status(const T& k) const {
        const auto p = no_doubles_.find(k);
        if (p == no_doubles_.end()) {
            return StatusType::UNQUEUED;
        }

        return p->second ? StatusType::PROCESSED : StatusType::QUEUED;
    }

    /// @brief Pop the top value in the queue.
    void pop() {
        pq_.pop();
    }

    /// @brief Top value in the queue.
    /// @return Top value
    std::pair<D, T> top() {
        return pq_.top();
    }

    /// @brief Bookkeeping for handling the next element.
    std::pair<D, T> next() {
        auto v = top();
        mark_as_processed(v.second);
        pop();
        return v;
    }

    /// @brief Empty and does bookkeeping. Call this before every next() call.
    bool empty() {
        sync();
        return pq_.empty();
    }

    /// @brief Add to be investigated
    /// @param d distance. Metric on which the priority queue is based on.
    ///        It returns always the minimum of the introduced values.
    /// @param id key value. Usually a tet_id.
    void push(const D& d, const T& id) {
        const auto stat = status(id);
        if (stat == StatusType::UNQUEUED) {
            pq_.emplace(d, id);
            no_doubles_[id] = false;
            update_min_unqueued_value();
        }
    }

    /// @brief Override the priority queue behavior and mark a key as already investigated.
    /// Previous status of the key in the queue does not matter. It will be treated as PROCESSED
    /// from now on.
    /// @param id Key value. Usually a tet_id
    void mark_as_processed(const T& id) {
        if (status(id) != StatusType::PROCESSED) {
            no_doubles_[id] = true;
            update_min_unqueued_value();
            ++size_processed_;
        }
    }

    /// @brief clear.
    void clear(const T& min_unqueued_value = T(0)) {
        min_unqueued_value_ = min_unqueued_value;
        pq_.clear();
        no_doubles_.clear();
        size_processed_ = 0;
    }

    /// @brief Get the smallest key that was not yet added to the priority queue.
    /// @return Key value. Usually a tet_id
    T get_min_unqueued_value() const {
        return min_unqueued_value_;
    }

    std::size_t size_processed() const {
        return size_processed_;
    }


  private:
    /// @brief Update the min_unqueued_value.
    void update_min_unqueued_value() {
        while (min_unqueued_value_ < T(int(no_doubles_.size())) &&
               status(min_unqueued_value_) != StatusType::UNQUEUED) {
            ++min_unqueued_value_;
        }
    }

    /// @brief Lazy sync of the priority queue the status of the keys. Keys that have been marked as
    /// PROCESSED are discarded without being processed.
    void sync() {
        while (!pq_.empty() && status(top().second) == StatusType::PROCESSED) {
            pq_.pop();
        }
    }

    /// @brief number of already processed keys.
    std::size_t size_processed_ = 0;
    /// @brief min unqueued value.
    T min_unqueued_value_;
    /// @brief priority queue.
    std::priority_queue<std::pair<D, T>, std::vector<std::pair<D, T>>, std::greater<>> pq_;
    /// @brief tagging system to keep track of the key statuses.
    std::unordered_map<T, bool> no_doubles_;
};


}  // namespace steps::util