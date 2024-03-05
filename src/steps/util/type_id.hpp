#pragma once

#include <cassert>
#include <limits>
#include <type_traits>
#include <utility>

namespace steps::util {

/**
 * \brief Data structure to store a pair (type, id) by
 * reusing the unused data space of id to store type.
 *
 * More specifically, it uses the \a EnumBits most significant bits
 * of Id to store Enum.
 *
 * It assumes that the identifier type is big enough to store
 * both all possible ids and types.
 *
 * For instance, if the values of an identifier are between 0 and 1 billion,
 * only 30 bits are required to encode them (2^31-1 = 1073741823).
 * The smallest appropriate integral type size is 4 bytes.
 * It means the 2 MSB can be reused for something else that can have
 * 4 values. Such type can be instantiated with:
 *
 *     enum class Type { A, B, C, D };
 *     using MyTypId = CompactType<Type, 2, unsigned>;
 *
 * \tparam Enum enum class representing the type
 * \tparam EnumBits number of bits required to encode all values of Enum
 * \tparam Id integer type of the identifier
 */
template <class Enum, int EnumBits, class Id>
class CompactTypeId {
  public:
    /// integer type of the identifier
    using id_t = Id;
    /// internal storage type
    using data_t = typename std::make_unsigned<Id>::type;

    static_assert(std::is_integral<id_t>::value, "Only integral types are supported for id_t");

    constexpr explicit CompactTypeId(data_t data)
        : data_(data) {}

    constexpr CompactTypeId(Enum type, id_t id)
        : data_(encode(type, id)) {
        assert(fits(type));
        assert(fits(id));
    }

    /**
     * \{
     * \name Accessors
     */

    constexpr Enum type() const noexcept {
        return static_cast<Enum>(data_ >> num_right_shifts_for_type());
    }

    constexpr id_t id() const noexcept {
        return static_cast<id_t>(data_ & id_mask());
    }

    constexpr data_t data() const noexcept {
        return data_;
    }

    /** \} */

    inline bool operator<(const CompactTypeId& rhs) const noexcept {
        return data_ < rhs.data_;
    }

    /**
     * \return an unsigned value with the same size than \a id_t
     * where \a type is written in MSB and \a id in LSB
     */
    static constexpr data_t encode(Enum type, id_t id) {
        return static_cast<data_t>(type) << num_right_shifts_for_type() | static_cast<data_t>(id);
    }

    /// \return the number of bits allocated to encode the type
    static constexpr int num_bits_for_type() {
        return EnumBits;
    }

    /// \return the total number of bits available to encode both type and id
    static constexpr int num_bits_for_data() {
        return std::numeric_limits<data_t>::digits;
    }

    /// ensure that type fits into the number of allocated bits
    static constexpr bool fits(Enum type) {
        const auto type_value = static_cast<data_t>(type);
        return (type_value & type_mask()) == type_value;
    }

    /// ensure that id fits into the number of allocated bits
    static constexpr bool fits(id_t id) {
        return (static_cast<data_t>(id) & id_mask()) == static_cast<data_t>(id);
    }

  private:
    static constexpr int num_right_shifts_for_type() {
        return num_bits_for_data() - num_bits_for_type();
    }

    static constexpr data_t type_mask() {
        return ~data_t{} >> num_right_shifts_for_type();
    }
    static constexpr data_t id_mask() {
        return ~data_t{} >> num_bits_for_type();
    }

    const data_t data_;
};

}  // namespace steps::util
