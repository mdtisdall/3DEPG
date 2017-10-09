#ifndef MKLAllocator_h
#define MKLAllocator_h

// Allocator that uses mkl_malloc/mkl_free.
//
// Modified from https://gist.github.com/donny-dont/1471329 
// 
// Modified from the Mallocator from Stephan T. Lavavej.
// <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>

#include <cstddef>
#include <stdexcept>
#include <new>

#include "mkl.h"

template <typename T >
class mkl_allocator
{
	public:
 
		// The following will be the same for virtually all allocators.
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T& reference;
		typedef const T& const_reference;
		typedef T value_type;
		typedef std::size_t size_type;
		typedef ptrdiff_t difference_type;
 
		T * address(T& r) const
		{
			return &r;
		}
 
		const T * address(const T& s) const
		{
			return &s;
		}
 
		std::size_t max_size() const
		{
			// The following has been carefully written to be independent of
			// the definition of size_t and to avoid signed/unsigned warnings.
			return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
		}
 
 
		// The following must be the same for all allocators.
		template <typename U>
		struct rebind
		{
			typedef mkl_allocator<U> other;
		} ;
 
		bool operator!=(const mkl_allocator& other) const
		{
			return !(*this == other);
		}
 
		void construct(T * const p, const T& t) const
		{
			void * const pv = static_cast<void *>(p);
 
			new (pv) T(t);
		}
 
		void destroy(T * const p) const
		{
			p->~T();
		}
 
		// Returns true if and only if storage allocated from *this
		// can be deallocated from other, and vice versa.
		// Always returns true for stateless allocators.
		bool operator==(const mkl_allocator&) const
		{
			return true;
		}
 
 
		// Default constructor, copy constructor, rebinding constructor, and destructor.
		// Empty for stateless allocators.
		mkl_allocator() { }
 
		mkl_allocator(const mkl_allocator&) { }
 
		template <typename U> mkl_allocator(const mkl_allocator<U>&) { }
 
		~mkl_allocator() { }
 
		void* allocate_inner(const std::size_t n) const {
      return mkl_malloc(n * sizeof(T), 64);
    }
 
		// The following will be different for each allocator.
		T * allocate(const std::size_t n) const
		{
			// The return value of allocate(0) is unspecified.
			// Mallocator returns NULL in order to avoid depending
			// on mkl_malloc(0)'s implementation-defined behavior
			// (the implementation can define malloc(0) to return NULL,
			// in which case the bad_alloc check below would fire).
			// All allocators can return NULL in this case.
			if (n == 0) {
				return NULL;
			}
 
			// All allocators should contain an integer overflow check.
			// The Standardization Committee recommends that std::length_error
			// be thrown in the case of integer overflow.
			if (n > max_size())
			{
				throw std::length_error("mkl_allocator<T>::allocate() - Integer overflow.");
			}
 
			// Mallocator wraps malloc().
			void * const pv = allocate_inner(n);
 
			// Allocators should throw std::bad_alloc in the case of memory allocation failure.
			if (pv == NULL)
			{
				throw std::bad_alloc();
			}
 
			return static_cast<T *>(pv);
		}


		void deallocate_inner(T * const p) const {
      mkl_free(p);
    }

		void deallocate(T * const p, const std::size_t) const
		{
			deallocate_inner(p);
		}
 
 
		// The following will be the same for all allocators that ignore hints.
		template <typename U>
		T * allocate(const std::size_t n, const U * /* const hint */) const
		{
			return allocate(n);
		}
 
 
		// Allocators are not required to be assignable, so
		// all allocators should have a private unimplemented
		// assignment operator. Note that this will trigger the
		// off-by-default (enabled under /Wall) warning C4626
		// "assignment operator could not be generated because a
		// base class assignment operator is inaccessible" within
		// the STL headers, but that warning is useless.
	private:
		mkl_allocator& operator=(const mkl_allocator&);
};

#endif
