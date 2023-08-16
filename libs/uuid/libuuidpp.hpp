//
//  libuuidcpp.hpp
//  C++ Class wrapper around libuuid
//
//  Boost Software License - Version 1.0 - August 17th, 2003
//
//  Permission is hereby granted, free of charge, to any person or organization
//  obtaining a copy of the software and accompanying documentation covered by
//  this license (the "Software") to use, reproduce, display, distribute,
//  execute, and transmit the Software, and to prepare derivative works of the
//  Software, and to permit third-parties to whom the Software is furnished to
//  do so, all subject to the following:
//
//  The copyright notices in the Software and this entire statement, including
//  the above license grant, this restriction and the following disclaimer,
//  must be included in all copies of the Software, in whole or in part, and
//  all derivative works of the Software, unless such copies or derivative
//  works are solely in the form of machine-executable object code generated by
//  a source language processor.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
//  FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
//  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//

#pragma once

#include <uuid/uuid.h>
#include <string>
#include <array>

#ifndef __APPLE__
#include <string.h>
/// 37 characters for a standard guid (including the zero terminator as per Darwin)
typedef char __uuid_string_t[37];
typedef __uuid_string_t  uuid_string_t;
#endif

// Compile-time assert if any expected types are not what we expect
static_assert(16 == sizeof(uuid_t), "uuid_t is of unexpected size");
static_assert(37 == sizeof(uuid_string_t), "Unexpected uuid_string_t size");

// MARK: - Typesafe enum bitfield support
//  - See https://www.justsoftwaresolutions.co.uk/cplusplus/using-enum-classes-as-bitfields.html

#include <type_traits>
#include <stdexcept>
#include <cstdint>

template<typename E>
struct enable_bitmask_operators{
	static const bool enable=false;
};

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable,E>::type
operator|(E lhs,E rhs){
	typedef typename std::underlying_type<E>::type underlying;
	return static_cast<E>(
		static_cast<underlying>(lhs) | static_cast<underlying>(rhs));
}

template<typename E>
typename std::enable_if<enable_bitmask_operators<E>::enable,E>::type
operator&(E lhs,E rhs){
	typedef typename std::underlying_type<E>::type underlying;
	return static_cast<E>(
		static_cast<underlying>(lhs) & static_cast<underlying>(rhs));
}

// MARK: - libuuidpp definition

namespace libuuidpp {

// Runtime exceptions
namespace exception {
	/// Exception for class constructor failure with bad (invalid) guids
	struct invalid_uuid: public std::runtime_error {
		explicit invalid_uuid(const char* guid)
		: std::runtime_error(std::string("libuuidpp exception: invalid guid '") + guid + "'")
		{}
	};
}

/// UUID class
class uuid {

public:

	/// Constant definition of the 'nil' uuid
	static const uuid nil;

	/// Binary representation class for the uuid
	typedef std::array<unsigned char, sizeof(uuid_t)> binary;

	// MARK: Static Creators

	/// Create a random uuid object.
	inline static uuid random() {
		return uuid(libuuidpp::uuid::Random);
	}

	// MARK: Validation

	/// Is the specific string a valid guid?
	inline static bool is_valid(const char* uuidString) {
		uuid_t tmp;
		return internal_create(uuidString, tmp);
	}

	/// Is the specific string a valid guid?
	inline static bool is_valid(const std::string& uuidString) {
		return is_valid(uuidString.c_str());
	}

	// MARK: Constructors

	/// Default constructor creates a nil guid
	inline uuid() {
		internal_set(libuuidpp::uuid::Nil);
	}

	/// Create a uuid object from a raw libuuid `uuid_t`
	inline explicit uuid(const uuid_t& rawValue) {
		::uuid_copy(_rawuuid, rawValue);
	}

	/// Create a uuid object from raw binary data
	inline explicit uuid(const binary& data) {
		set(data);
	}

	/**
	 Create a new uuid object using the provided guid string

	 @param guidString the string version the guid
	 @throws uuid::invalid if the guid is invalid
	 */
	inline explicit uuid(const char* guidString) {
		if (!set(guidString)) {
			throw exception::invalid_uuid((guidString != nullptr) ? guidString : "NULL");
		}
	}

	/// Convenience constructor for std::string
	inline explicit uuid(const std::string& guidString)
		: uuid(guidString.c_str()) {
	}

	/// Copy constructor
	inline uuid(const uuid& right) {
		::uuid_copy(_rawuuid, right._rawuuid);
	}

	// MARK: Setters

	/**
	 Set a uuid using the provided guid string

	 @param uuidString the string version the guid
	 @returns true if the guid string is valid and assignable, false otherwise
	 */
	inline bool set(const char* uuidString) {
		return internal_create(uuidString, _rawuuid);
	}

	inline bool set(const std::string& guidString) {
		return set(guidString.c_str());
	}

	inline void set(const binary& data) {
		memcpy(_rawuuid, data.data(), data.size());
	}

	/// Set a random guid
	inline void set_random() {
		::uuid_generate_random(_rawuuid);
	}
	/// Set a nil guid
	inline void set_nil() {
		::uuid_clear(_rawuuid);
	}

	// MARK: Hash functions

	/**
	 Return a hash of the uuid

	 Generate a hash value for the uuid for use in unordered_set and unordered_map.

	 Uses Fowler–Noll–Vo FNV-1a hash function (https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function)

	 @return the hash value
	 */
	[[nodiscard]] std::size_t hash() const;

	// MARK: Comparators

	[[nodiscard]] inline bool is_nil() const {
		return ::uuid_is_null(_rawuuid);
	}
	inline bool operator!=(const uuid& right) const {
		return ::uuid_compare(_rawuuid, right._rawuuid) != 0;
	}
	inline bool operator==(const uuid& right) const {
		return ::uuid_compare(_rawuuid, right._rawuuid) == 0;
	}
	inline bool operator<(const uuid& right) const {
		return ::uuid_compare(_rawuuid, right._rawuuid) < 0;
	}
	inline bool operator>(const uuid& right) const {
		return ::uuid_compare(_rawuuid, right._rawuuid) > 0;
	}

	/// The available formattings for string output (bitfield)
	enum class formatting {
		/// Standard formatting (uppercased, no brackets)
		standard = 0,
		/// Lowercase output
		lowercase = 1 << 0,
		/// Use Microsoft format for output (leading and trailing bracket)
		brackets = 1 << 1
	};

	/**
	 Returns a string representation of the uuid

	 @param format The format type, lowercase, uppercase, bracketted etc.
	 @return String representation of the uuid
	 */
	[[nodiscard]] std::string string(formatting format = formatting::standard) const;

	/// Copy out the raw uuid_t value
	[[nodiscard]] inline uuid::binary data() const {
		uuid::binary tmp;
		::memcpy(tmp.data(), _rawuuid, sizeof(uuid_t));
		return tmp;
	}

// MARK: - Private implementation

private:

	/// Internal creation type
	typedef enum CreationState {
		/// Fill the created guid with a nil guid
		Nil = 0,
		/// Assign a random value to the created guid
		Random = 1
	} CreationState;

	// EB758F6F-E2B0-46F1-BCDF-4162A6EB11D9
	// Note that uuid_string_t INCLUDES the zero terminator
	static const size_t STR_GUID_LEN = sizeof(uuid_string_t);

	// {EB758F6F-E2B0-46F1-BCDF-4162A6EB11D9}
	// Microsoft generated guids can have brackets around them.
	// NOTE: The length here is only 1 more (instead of the expected 2) because the
	//       defined size of `uuid_string_t` includes the NULL terminator.
	static const size_t MS_GUID_LEN = STR_GUID_LEN + 1;

	inline static uuid internal_create(uuid::CreationState creation = uuid::Nil) {
		return uuid(creation);
	}

	inline static bool is_microsoft_formatted(const char* guid) {
		return ((::strnlen(guid, MS_GUID_LEN + 2) == MS_GUID_LEN)
				&& guid[0] == '{'
				&& guid[STR_GUID_LEN] == '}');
	}

	static bool internal_microsoft_create(const char* uuidString, uuid_t& result);

	inline static bool internal_uuid_create(const char* uuidString, uuid_t& result) {
		return (::uuid_parse(uuidString, result) == 0);
	}

	static bool internal_create(const char* uuidString, uuid_t& result);

	inline void internal_set(CreationState creation = uuid::Nil) {
		(creation == uuid::Random) ? set_random() : set_nil();
	}

	inline explicit uuid(CreationState creation) {
		internal_set(creation);
	}

	/// The raw underlying libuuid data
	::uuid_t _rawuuid{};
};
};

/// Install the hashing function for libuuidpp::uuid into the stdlib
namespace std {
	template<> struct hash<libuuidpp::uuid> {
		std::size_t operator()(libuuidpp::uuid const& s) const noexcept {
			return s.hash();
		}
	};
}

template<>
struct enable_bitmask_operators<libuuidpp::uuid::formatting> {
	static const bool enable=true;
};

// MARK: - libuuidpp Implementations

namespace libuuidpp {

const uuid uuid::nil = uuid::internal_create(libuuidpp::uuid::Nil);

std::string uuid::string(libuuidpp::uuid::formatting format) const {

	std::string result;
	result.reserve(MS_GUID_LEN);

	uuid_string_t strVal;

	bool isBracketed = (format & formatting::brackets) == formatting::brackets;
	bool isLowercased = (format & formatting::lowercase) == formatting::lowercase;

	if (isBracketed) {
		result += "{";
	}

	if (isLowercased) {
		::uuid_unparse_lower(_rawuuid, strVal);
	}
	else {
		::uuid_unparse_upper(_rawuuid, strVal);
	}
	result += std::string(strVal);

	if (isBracketed) {
		result += "}";
	}

	return result;
}

std::size_t uuid::hash() const {

#if INTPTR_MAX == INT64_MAX		// 64-bit
	static const std::uint64_t FNV_offset_basis = 14695981039346656037UL;
	static const std::uint64_t FNV_prime = 1099511628211;
#elif INTPTR_MAX == INT32_MAX	// 32-bit
	static const std::uint32_t FNV_offset_basis = 2166136261;
	static const std::uint32_t FNV_prime = 16777619;
#else
#error Not a supported architecture (32bit or 64bit)
#endif

	std::size_t result = FNV_offset_basis;
	for (std::size_t offset = 0; offset < sizeof(uuid_t); offset++) {
		result ^= _rawuuid[offset];
		result *= FNV_prime;
	}
	return result;
}

bool uuid::internal_create(const char* uuidString, uuid_t& result) {
	if (uuidString != nullptr) {
		if (internal_microsoft_create(uuidString, result)) {
			return true;
		}
		else if (internal_uuid_create(uuidString, result)) {
			return true;
		}
	}
	::uuid_clear(result);
	return false;
}

bool uuid::internal_microsoft_create(const char* uuidString, uuid_t& result) {
	if (is_microsoft_formatted(uuidString)) {
		char temp[STR_GUID_LEN];
		::memset(temp, 0, STR_GUID_LEN);
		::memcpy(temp, uuidString+1, STR_GUID_LEN - 1);
		if (::uuid_parse(temp, result) != 0) {
			return false;
		}
		return true;
	}
	return false;
}

}
