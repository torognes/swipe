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

// libFuzzer harness of the ASN.1 header parser (asnparse.cc): each
// input is the header of one database sequence (a Blast-def-line-set,
// as makeblastdb writes it in the .phr/.nhr files), handed to the
// three entry points that swipe calls for each hit or sequence.
//
// Build and run with fuzz/run_fuzzer.sh (clang only).
//
// fatal() exits in swipe; here it jumps back to the harness instead,
// so that one rejected header does not end the fuzzing run. The jump
// skips the destructors of the parser's local strings and vectors:
// the harness keeps a list of the blocks allocated by operator new
// while an input is parsed, and frees those still allocated once the
// parser is gone (without that, the leaks of millions of rejected
// headers exhausted the memory and libFuzzer reported out-of-memory
// errors, blamed on whatever input was running). Each input gets its
// own parser, so that a crash reproduces from its input alone, and
// LeakSanitizer stays on: a leak on a path without fatal() is real.

#include "../swipe.h"
#include <climits>  // LONG_MAX
#include <csetjmp>
#include <cstddef>  // std::size_t, std::max_align_t
#include <cstdint>  // std::uint8_t
#include <cstdio>  // std::fopen, std::fputc
#include <cstdlib>  // std::malloc, std::free
#include <new>  // std::bad_alloc
#include <string>

// the definitions that swipe.cc and hits.cc provide to asnparse.cc
FILE * out = nullptr;

namespace {

std::jmp_buf on_fatal;

// --- the blocks allocated while an input is parsed ---------------------------

// every block of operator new starts with a node (padded to keep the
// alignment of malloc); the nodes of the blocks allocated while an
// input is parsed form a doubly-linked list
struct alignas(std::max_align_t) Node
{
  Node * previous;
  Node * next;
  bool listed;
};

Node listed_blocks {&listed_blocks, &listed_blocks, false};
bool listing = false;

auto unlink(Node & node) -> void
{
  node.previous->next = node.next;
  node.next->previous = node.previous;
  node.listed = false;
}

// a block of operator new, after its node (nullptr if out of memory)
auto allocate(std::size_t const size) noexcept -> void *
{
  auto * const node = static_cast<Node *>(std::malloc(sizeof(Node) + size));
  if (node == nullptr)
  {
    return nullptr;
  }
  node->listed = listing;
  if (listing)
  {
    node->previous = &listed_blocks;
    node->next = listed_blocks.next;
    listed_blocks.next->previous = node;
    listed_blocks.next = node;
  }
  return node + 1;
}

// free the blocks still allocated (skipped destructors)
auto free_listed_blocks() -> void
{
  while (listed_blocks.next != &listed_blocks)
  {
    auto * const node = listed_blocks.next;
    unlink(*node);
    std::free(node);
  }
}

// every taxid passes (swipe without -x), or one in two (with -x)
auto accept_all(long /* taxid */) -> long
{
  return 1;
}

auto accept_odd(long const taxid) -> long
{
  return taxid % 2;
}

}  // namespace

auto operator new(std::size_t const size) -> void *
{
  auto * const block = allocate(size);
  if (block == nullptr)
  {
    throw std::bad_alloc();
  }
  return block;
}

auto operator delete(void * const block) noexcept -> void
{
  if (block == nullptr)
  {
    return;
  }
  auto * const node = static_cast<Node *>(block) - 1;
  if (node->listed)
  {
    unlink(*node);
  }
  std::free(node);
}

auto operator new[](std::size_t const size) -> void *
{
  return operator new(size);
}

auto operator delete[](void * const block) noexcept -> void
{
  operator delete(block);
}

// the other forms called by libFuzzer: nothrow, and sized (C++14,
// -fsized-deallocation). The aligned forms (C++17) are left to the
// sanitizer runtime: an aligned new is always paired with an aligned
// delete.
auto operator new(std::size_t const size, std::nothrow_t const & /* tag */) noexcept -> void *
{
  return allocate(size);
}

auto operator new[](std::size_t const size, std::nothrow_t const & /* tag */) noexcept -> void *
{
  return allocate(size);
}

auto operator delete(void * const block, std::nothrow_t const & /* tag */) noexcept -> void
{
  operator delete(block);
}

auto operator delete[](void * const block, std::nothrow_t const & /* tag */) noexcept -> void
{
  operator delete(block);
}

auto operator delete(void * const block, std::size_t /* size */) noexcept -> void
{
  operator delete(block);
}

auto operator delete[](void * const block, std::size_t /* size */) noexcept -> void
{
  operator delete(block);
}

[[noreturn]] auto fatal(char const * /* message */) noexcept -> void
{
  std::longjmp(on_fatal, 1);
}

[[noreturn]] auto fatal(std::string const & /* message */) noexcept -> void
{
  std::longjmp(on_fatal, 1);
}

auto xml_putc(char const symbol) noexcept -> void
{
  std::fputc(symbol, out);
}

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const * data,
                                       std::size_t size) -> int;

extern "C" auto LLVMFuzzerTestOneInput(std::uint8_t const * data,
                                       std::size_t const size) -> int
{
  if (out == nullptr)
  {
    out = std::fopen("/dev/null", "w");
  }

  // the layouts of the output formats (-m 0, -m 8, -m 9, -m 99)
  HeaderLayout plain;
  HeaderLayout wrapped;
  wrapped.indent = 2;
  wrapped.maxlen = 40;
  wrapped.linelen = 30;
  wrapped.maxdeflines = 3;
  HeaderLayout identifier;
  identifier.text = DeflineText::identifier;
  identifier.escaping = Escaping::xml;
  identifier.show_gis = 1;

  listing = true;
  {
    auto const parser = parser_create(1);  // with -T (show taxids)

    // a copy of exactly size bytes, so that AddressSanitizer reports
    // any read past the end of the header
    std::string const copy(reinterpret_cast<char const *>(data), size);
    auto const header = View<char>(copy.data(), copy.size());

    // volatile: modified between setjmp() and longjmp()
    for (volatile int step = 0; step < 6; step = step + 1)
    {
      if (setjmp(on_fatal) != 0)
      {
        continue;  // header rejected by fatal()
      }
      switch (step)
      {
      case 0:
        static_cast<void>(parse_getdeflines(*parser, header, 0, accept_all, 0));
        break;
      case 1:
        static_cast<void>(parse_getdeflines(*parser, header, 1, accept_odd, 1));
        break;
      case 2:
        static_cast<void>(parse_getdeflinecount(*parser, header, 0, accept_odd));
        break;
      case 3:
        parse_header(*parser, header, 0, accept_all, plain);
        break;
      case 4:
        parse_header(*parser, header, 0, accept_all, wrapped);
        break;
      default:
        parse_header(*parser, header, 0, accept_all, identifier);
        break;
      }
    }
  }
  listing = false;
  free_listed_blocks();
  return 0;
}
