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

// imported from vsearch (src/utils/fatal_allocator.hpp, commit
// a7cb99e6), with swipe's xmalloc() and free()

#ifndef SWIPE_FATAL_ALLOCATOR_H
#define SWIPE_FATAL_ALLOCATOR_H

#include <cstddef>  // std::size_t
#include <cstdlib>  // std::free
#include <vector>


// a 16-byte aligned block, or fatal() when out of memory (swipe.cc)
auto xmalloc(std::size_t size) -> void *;


/* A minimal standard-library allocator that obtains memory through xmalloc and
   releases it through free. xmalloc reports an out-of-memory condition through
   fatal() (a clean message, then exit) and never returns null, so an exhausted
   allocation ends the program exactly as the raw buffers always did, instead
   of throwing std::bad_alloc -- which, under -fno-exceptions, would
   std::terminate. swipe's xmalloc also aligns every block on 16 bytes, as the
   SIMD buffers need. Stateless: any two instances compare equal, so a
   container may move or swap its storage freely. C++11 std::allocator_traits
   synthesises the rest of the allocator interface (rebind, construct,
   max_size, ...). */
template <class Element>
struct FatalAllocator
{
  using value_type = Element;

  FatalAllocator() = default;

  // enables a container to rebind the allocator for its internal element type
  template <class Other>
  FatalAllocator(FatalAllocator<Other> const & /*other*/) noexcept {}  // NOLINT: implicit by design

  auto allocate(std::size_t const count) -> Element *
  {
    // a container never requests more than max_size() elements, so the product
    // cannot overflow std::size_t; xmalloc fatals rather than returning null
    return static_cast<Element *>(xmalloc(count * sizeof(Element)));
  }

  auto deallocate(Element * const pointer, std::size_t const /*count*/) noexcept -> void
  {
    std::free(static_cast<void *>(pointer));
  }
};


template <class Left, class Right>
auto operator==(FatalAllocator<Left> const & /*lhs*/, FatalAllocator<Right> const & /*rhs*/) noexcept -> bool
{
  return true;
}

template <class Left, class Right>
auto operator!=(FatalAllocator<Left> const & /*lhs*/, FatalAllocator<Right> const & /*rhs*/) noexcept -> bool
{
  return false;
}

// a std::vector whose storage comes from xmalloc (16-byte aligned; out of
// memory: fatal)
template <class Element>
using Buffer = std::vector<Element, FatalAllocator<Element>>;

#endif  // SWIPE_FATAL_ALLOCATOR_H
