/*
 * ============================================================================
 *
 *       Filename:  cqf.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017-10-26 11:50:04 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (), ppandey@cs.stonybrook.edu
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#ifndef _CQF_H_
#define _CQF_H_

#include <iostream>
#include <cassert>
#include <unordered_set>

#include <inttypes.h>
#include <string.h>
#include <math.h>

#include "gqf/gqf.h"
#include "gqf/gqf_file.h"
#include "util.h"

#define NUM_HASH_BITS 24
#define NUM_Q_BITS 16

enum readmode {
	MMAP,
	FREAD
};

template <class key_obj>
class CQF {
	public:
		CQF();
		CQF(uint64_t q_bits, uint64_t key_bits, enum lockingmode lock, enum
				hashmode hash, uint32_t seed);
		CQF(std::string& filename, enum lockingmode lock, enum readmode flag);
		CQF(const CQF<key_obj>& copy_cqf);

		bool insert(const key_obj& k);

		/* Will return the count. */
		uint64_t query(const key_obj& k);

		uint64_t inner_prod(const CQF<key_obj>& in_cqf);

		void serialize(std::string filename) {
			qf_serialize(&cqf, filename.c_str());
		}

		void set_auto_resize(void) { qf_set_auto_resize(&cqf); }

		bool is_exact(void);

		const QF* get_cqf(void) const { return &cqf; }
		uint64_t range(void) const { return cqf.metadata->range; }
		uint32_t seed(void) const { return cqf.metadata->seed; }
		uint32_t keybits(void) const { return cqf.metadata->key_bits; }
		uint64_t total_elts(void) const { return cqf.metadata->nelts; }
		uint64_t dist_elts(void) const { return cqf.metadata->ndistinct_elts; }
		//uint64_t set_size(void) const { return set.size(); }
		void reset(void) { qf_reset(&cqf); }

		void dump_metadata(void) const { DEBUG_DUMP(&cqf); }

		void drop_pages(uint64_t cur);

		class Iterator {
			public:
				Iterator(QFi it);
				key_obj operator*(void) const;
				void operator++(void);
				bool done(void) const;

				QFi iter;
			private:
				uint64_t end_hash;
		};

		Iterator begin(void) const;
		Iterator end(void) const;

	private:
		QF cqf;
		//std::unordered_set<uint64_t> set;
};

class KeyObject {
	public:
		KeyObject() : key(0), value(0), count(0) {};

		KeyObject(uint64_t k, uint64_t v, uint64_t c) : key(k),
		value(v), count(c) {};

		KeyObject(const KeyObject& k) : key(k.key), value(k.value), count(k.count) {};

		bool operator==(KeyObject k) { return key == k.key; }

		uint64_t key;
		uint64_t value;
		uint64_t count;
};

template <class key_obj>
CQF<key_obj>::CQF() {
	if (!qf_malloc(&cqf, 1ULL << NUM_Q_BITS, NUM_HASH_BITS, 0, LOCKS_REQUIRED,
								 DEFAULT, 2038074761)) {
		ERROR("Can't allocate the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(uint64_t q_bits, uint64_t key_bits, enum lockingmode lock,
									enum hashmode hash, uint32_t seed) {
	if (!qf_malloc(&cqf, 1ULL << q_bits, key_bits, 0, lock, hash, 2038074761)) {
		ERROR("Can't allocate the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(std::string& filename, enum lockingmode lock, enum readmode
									flag) {
	uint64_t size = 0;
	if (flag == MMAP)
	 size = qf_usefile(&cqf, lock, filename.c_str());
	else
		size = qf_deserialize(&cqf, lock, filename.c_str());

	if (size == 0) {
		ERROR("Can't read/deserialize the CQF");
		exit(EXIT_FAILURE);
	}
}

template <class key_obj>
CQF<key_obj>::CQF(const CQF<key_obj>& copy_cqf) {
	memcpy(&cqf, copy_cqf.get_cqf(), sizeof(QF));
}

template <class key_obj>
bool CQF<key_obj>::insert(const key_obj& k) {
	return qf_insert(&cqf, k.key, k.value, k.count);
	// To validate the CQF
	//set.insert(k.key);
}

template <class key_obj>
uint64_t CQF<key_obj>::query(const key_obj& k) {
	return qf_count_key_value(&cqf, k.key, k.value);
}

template <class key_obj>
uint64_t CQF<key_obj>::inner_prod(const CQF<key_obj>& in_cqf) {
	return qf_inner_product(&cqf, in_cqf.get_cqf());
}

template <class key_obj>
bool CQF<key_obj>::is_exact(void) {
	if (cqf.metadata->hash_mode == INVERTIBLE)
		return true;
	return false;
}

template <class key_obj>
CQF<key_obj>::Iterator::Iterator(QFi it)
	: iter(it) {};

template <class key_obj>
key_obj CQF<key_obj>::Iterator::operator*(void) const {
	uint64_t key = 0, value = 0, count = 0;
	qfi_get(&iter, &key, &value, &count);
	return key_obj(key, value, count);
}

template<class key_obj>
void CQF<key_obj>::Iterator::operator++(void) {
	qfi_next(&iter);
}

/* Currently, the iterator only traverses forward. So, we only need to check
 * the right side limit.
 */
template<class key_obj>
bool CQF<key_obj>::Iterator::done(void) const {
	return qfi_end(&iter);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::begin(void) const {
	QFi qfi;
	qf_iterator(&this->cqf, &qfi, 0);
	return Iterator(qfi);
}

template<class key_obj>
typename CQF<key_obj>::Iterator CQF<key_obj>::end(void) const {
	QFi qfi;
	qf_iterator(&this->cqf, &qfi, 0xffffffffffffffff);
	return Iterator(qfi, UINT64_MAX);
}

#endif
