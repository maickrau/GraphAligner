#ifndef ArenaAllocator_h
#define ArenaAllocator_h

#include <memory>

#ifdef MEMORYVERBOSE
#include <iostream>
#endif


#ifdef MEMORYASSERTS
#include "ThreadReadAssertion.h"
#define memassert(expression) assert(expression)
#else
#define memassert(ignore) ((void)0)
#endif


namespace ArenaAllocator 
{
		static constexpr size_t BLOCKSIZE = 16*1024*1024; //arbitrarily 16 Mb

	struct MemBlock
	{
	public:
#ifdef MEMORYVERBOSE
		~MemBlock()
		{
			std::cerr << "deallocate block " << std::endl;
			std::cerr << "allocates " << numAllocates << std::endl;
			std::cerr << "bytes " << byteAllocates << std::endl;
		}
		size_t byteAllocates;
		size_t numAllocates;
#endif
		std::vector<char> start;
		size_t sizeleft;
		char* pos;
		std::unique_ptr<MemBlock> previous;
	};

	class MemChain
	{
	public:
		std::unique_ptr<MemBlock> head;
	};

	std::unique_ptr<MemBlock> allocateBlock()
	{
		std::unique_ptr<MemBlock> result = std::make_unique<MemBlock>();
		result->start.resize(BLOCKSIZE);
#ifdef MEMORYVERBOSE
		std::cerr << "allocate block" << std::endl;
#endif
		// memset((void*)result->start, 0, BLOCKSIZE);
		result->pos = result->start.data();
		result->sizeleft = result->start.size();
		result->previous = nullptr;
#ifdef MEMORYVERBOSE
		result->numAllocates = 0;
		result->byteAllocates = 0;
#endif
		return result;
	}

	std::shared_ptr<MemChain> makeChain()
	{
#ifdef MEMORYVERBOSE
		std::cerr << "allocate chain" << std::endl;
#endif
		auto result = std::make_shared<MemChain>();
		result->head = allocateBlock();
		return result;
	}

	thread_local std::shared_ptr<MemChain> sharedChain = makeChain();

	void ResetAllocator()
	{
		sharedChain = makeChain();
	}

	template <typename T>
	class ArenaAllocator
	{
	public:
		using allocator_type = T;
		using value_type = T;
		ArenaAllocator(const ArenaAllocator& other) :
		currentChain(other.currentChain)
		{
		}
		ArenaAllocator& operator=(const ArenaAllocator& other)
		{
			if (&other == this) return *this;
			currentChain = other.currentChain;
			return *this;
		}
		template <typename U>
		ArenaAllocator(const ArenaAllocator<U>& other) :
		currentChain(other.currentChain)
		{
		}
		template <typename U>
		ArenaAllocator& operator=(const ArenaAllocator<U>& other)
		{
			currentChain = other.currentChain;
			return *this;
		}
		bool operator==(const ArenaAllocator& other) const noexcept
		{
			if (currentChain == nullptr || other.currentChain == nullptr) return false;
			return currentChain == other.currentChain;
		}
		bool operator!=(const ArenaAllocator& other) const noexcept
		{
			return !(other == *this);
		}
		ArenaAllocator() :
		currentChain(sharedChain)
		{
		}
		~ArenaAllocator()
		{
		}
		T* address(T& item) const noexcept
		{
			return &item;
		}
		T* address(const T& item) const noexcept
		{
			return &item;
		}
		T* allocate(size_t count)
		{
			size_t size = count * sizeof(T);
			memassert(currentChain != nullptr);
			memassert(currentChain->head != nullptr);
			memassert(size < BLOCKSIZE);
			if (size >= currentChain->head->sizeleft)
			{
				std::unique_ptr<MemBlock> newBlock = allocateBlock();
				newBlock->previous = std::move(currentChain->head);
				currentChain->head = std::move(newBlock);
			}
#ifdef MEMORYVERBOSE
			currentChain->head->byteAllocates += size;
			currentChain->head->numAllocates += 1;
#endif
			T* result = (T*)currentChain->head->pos;
			memassert(currentChain->head->end - currentChain->head->pos >= size);
			memassert(currentChain->head->pos < currentChain->head->end);
			currentChain->head->pos += size;
			currentChain->head->sizeleft -= size;
			return result;
		}
		void deallocate(T* ptr, size_t count)
		{
		}
		size_t max_size() const noexcept
		{
			return BLOCKSIZE / sizeof(T);
		}
		//http://www.cplusplus.com/reference/memory/allocator/construct/
		template <class U, class... Args>
		void construct (U* p, Args&&... args)
		{
			::new ((void*)p) U (std::forward<Args>(args)...);
		}
		template <class U>
		void destroy (U* p)
		{
			p->~U();
		}
	private:
		std::shared_ptr<MemChain> currentChain;
		template <typename U>
		friend class ArenaAllocator;
	};

}

#endif
