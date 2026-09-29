/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2021 Torbjorn Rognes, University of Oslo,
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

// imported from vsearch (src/utils/element_order.hpp, commit 63acedba)

#ifndef SWIPE_ELEMENT_ORDER_H
#define SWIPE_ELEMENT_ORDER_H

#include <limits>  // std::numeric_limits
#include <type_traits>  // std::remove_cv


// How View and Span order their elements.
//
// The primary template defers to the element type's own operator<, which is
// what a bare std::lexicographical_compare would have done anyway.
//
// The specialization for char is the reason this trait exists: it orders bytes
// as unsigned char. That is what std::strcmp does, and what std::string does
// through std::char_traits<char>::lt -- but it is *not* what comparing char
// with operator< does, because char is signed on x86-64 and on the Windows
// target while ARM and PowerPC Linux default it to unsigned. Without the
// specialization, View<char>{} < View<char>{} would order any byte with its
// high bit set differently from std::string, differently from std::strcmp, and
// differently from one architecture to the next.
//
// That divergence is reachable and user-visible in vsearch: sequence headers
// are arbitrary bytes and routinely carry UTF-8 or Latin-1 (an accented author
// name, a non-ASCII locality), and the header is the tie-break of every
// --sortby* and --derep* ordering. Measured over 2692881 pairs from a byte
// alphabet straddling 0x80, the unspecialized comparison disagrees with
// std::strcmp on 1378400 of them where char is signed, and on none where it is
// unsigned.
//
// Only char is specialized. signed char and unsigned char say what they mean
// and are left to the primary template, the same line std::char_traits draws.

template <typename Type>
struct element_order_unqualified {
  static constexpr auto less(Type const & lhs, Type const & rhs) -> bool {
    return lhs < rhs;
  }

  // Three-way, with the sign convention of std::strcmp: negative if lhs sorts
  // first, positive if rhs does, zero if the two are equivalent. Expressed with
  // two less() calls because the primary template can assume nothing beyond a
  // strict weak ordering.
  static auto compare(Type const & lhs, Type const & rhs) -> int {
    if (less(lhs, rhs)) { return -1; }
    if (less(rhs, lhs)) { return +1; }
    return 0;
  }
};

// Both members are noexcept, unlike the primary template's: a conversion to
// unsigned char and an integral comparison are the whole computation, and
// neither can throw. This is the leaf that View<char>'s and Span<char>'s
// noexcept comparison members bottom out in (see span.hpp), so the guarantee
// is honest all the way down for the one specialization that matters.
template <>
struct element_order_unqualified<char> {
  static constexpr auto less(char const lhs, char const rhs) noexcept -> bool {
    return static_cast<unsigned char>(lhs) < static_cast<unsigned char>(rhs);
  }

  // One byte comparison rather than the primary template's two: char is
  // totally ordered once it is read as unsigned char.
  static auto compare(char const lhs, char const rhs) noexcept -> int {
    auto const lhs_byte = static_cast<unsigned char>(lhs);
    auto const rhs_byte = static_cast<unsigned char>(rhs);
    if (lhs_byte == rhs_byte) { return 0; }
    return (lhs_byte < rhs_byte) ? -1 : +1;
  }
};

// The entry point Span and View actually use, keyed on the cv-unqualified
// element type. The stripping is what keeps Span<char const> and
// View<char const> -- spellings both class templates accept, see comparable
// in span.hpp -- ordering their bytes as unsigned char. Without it those
// instantiations miss the specialization above and fall through to the primary
// template, which is the very divergence from std::strcmp that the
// specialization exists to prevent: 'A' < '\xc3' comes out false where char is
// signed, true where it is unsigned, and true for Span<char> on both.
//
// An alias rather than a class template deriving from the implementation: it
// adds no base class, and it makes the wrong extension point a compile error,
// since a future specialization cannot be written on an alias template and has
// to go on element_order_unqualified where the stripping still reaches it.
template <typename Type>
using element_order =
  element_order_unqualified<typename std::remove_cv<Type>::type>;

// Guards for the paragraph above; constexpr less() makes them free. Wrapped in
// a struct so that the byte constant does not reach the global namespace of
// every TU that includes this header, for the reason spelled out for
// max_ptrdiff in span.hpp.
struct element_order_guards {
  // The byte 'A' is ordered against, spelled as unsigned char's maximum rather
  // than as a literal so that it means "high bit set" on every target: it is
  // negative where char is signed (x86-64, Windows) and positive where it is
  // unsigned (ARM and PowerPC Linux), and both assertions hold either way.
  static constexpr char high_bit_byte =
    static_cast<char>(std::numeric_limits<unsigned char>::max());
  static_assert(element_order<char>::less('A', high_bit_byte),
                "Span<char> and View<char> must order bytes as unsigned char");
  static_assert(element_order<char const>::less('A', high_bit_byte),
                "a cv-qualified element type must order like a bare one");
};


// Functor wrapper, for the std algorithms that take a comparison object.
//
// No stripping of its own: it routes through element_order above, so
// element_less<char const> already compares like element_less<char>.
template <typename Type>
struct element_less {
  auto operator()(Type const & lhs, Type const & rhs) const -> bool {
    return element_order<Type>::less(lhs, rhs);
  }
};

#endif  // SWIPE_ELEMENT_ORDER_H
