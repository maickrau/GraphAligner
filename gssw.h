/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 01/30/13.
 *
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <smmintrin.h>

/*!	@typedef	structure of the query profile	*/
struct gssw_profile;
typedef struct gssw_profile gssw_profile;

// This file makes extensive use of SSE intrinsics, which operate on 128-bit
// values. These values are stored in the 128-bit special registers XMM0 through
// XMM7, as assigned by the compiler. In C, these values have type "__m128i",
// with the "i" meaning integer (because we do integer math on them). We have
// code to use these in alignments, either to store 16 single-byte values or 8
// double-byte values, depending on how large we're allowing the scores to get.

// We use j for indexes in the read (or "query"), and i for indexes in the
// reference (or "database"). Conceptually, the database runs rightwards along
// the top of the DP matrix, and the read runs downwards along the left.

// Note that this usage of i and j is BACKWARDS from the paper!

// Within each node, we perform a standard-ish Smith-Waterman alignment with
// affine gaps. This involves three matrices. H[i, j] gives the score of the
// very best alignment using read characters through j and reference characters
// through i, no matter what it ends in (match, mismatch, or gap). We also have
// E[i, j] which gives the score of the best such alignment ending in a gap in
// the read, and F[i, j], which give sthe score of the best such alignment
// ending in a gap in the reference.

// During alignment, we work on one column (i.e. reference character) of all the
// DP matrices at a time. This column is *not* divided top to bottom into SSE
// vectors. Instead, the column is broken up across SSE vectors in a "striped"
// format, as follows:

// Take the read, and break it into 16 (or 8 for word-scored alignment) equal-
// length pieces, padding if necessary. Line those pieces up left to right. Then
// run through *that* matrix in rows. Each row is an SSE 128-bit chunk, and is
// processed all at once.

// Another way to think of the layout is to have 16 or 8 cursors moving down the
// column in parallel. Each cursor has easy access to the contents of the cell
// above it, except for the first one (which is fixed with some fix-up logic
// called the "lazy F loop").

// Each node has a notion of a "seed", which carries the E information (about
// gaps in the read that can continue through the node boundary) and the H
// information (about the best alignments overall before we entered this node's
// matrix). F information is not carried over, because gaps in the reference
// can't span multiple reference nodes.
typedef struct {
    __m128i* pvE;
    __m128i* pvHStore;
} gssw_seed;


/*!	@typedef	structure of the alignment result
	@field	score1	the best alignment score
	@field	score2	sub-optimal alignment score
	@field	ref_begin1	0-based best alignment beginning position on reference;	ref_begin1 = -1 when the best alignment beginning
						position is not available
	@field	ref_end1	0-based best alignment ending position on reference
	@field	read_begin1	0-based best alignment beginning position on read; read_begin1 = -1 when the best alignment beginning
						position is not available
	@field	read_end1	0-based best alignment ending position on read
	@field	read_end2	0-based sub-optimal alignment ending position on read
*/
typedef struct {
	uint16_t score1;
	uint16_t score2;
	int32_t ref_begin1;
	int32_t ref_end1;
	int32_t	read_begin1;
	int32_t read_end1;
	int32_t ref_end2;
	// We need the seed information about H and E, which we get from the possible
	// previous nodes.
    gssw_seed seed;
    // This stores 1 if we did the alignment with byte scores, in which case the
    // matrix pointers below should point to an matrices composed of bytes. If
    // it is 0, the matrices are composed of 2-byte values.
    uint8_t is_byte;
    
    // We keep de-striped versions of all the alignment matrices, to allow
    // traceback.
    
    // Pointer to the H matrix for the alignment, which holds, at each position,
    // the score of the best alignment of the corresponding substrings, no
    // matter what it ends with.
    void* mH;
    // Pointer to the E matrix for the alignment, which holds scores for
    // alignments ending with gaps in the read.
    void* mE;
    // Pointer to the F matrix, which holds scores for alignments ending with
    // gaps in the reference.
    void* mF;
} gssw_align;

typedef struct {
	uint16_t score;
	int32_t ref;	 //0-based position
	int32_t read;    //alignment ending position on read, 0-based
} gssw_alignment_end;

// A CIGAR element is a number of repetitions of a CIGAR operation.
typedef struct {
    char type;
    uint32_t length;
} gssw_cigar_element;

// A CIGAR is an array of elements.
typedef struct {
    int32_t length;
    gssw_cigar_element* elements;
} gssw_cigar;

// A profile is a sort of exploded score matrix. We fill a whole read by
// reference matrix with the score you would get by matching a particular
// character in the read with a particular character in the reference. There's
// no data dependence between any entries, so we can just fill this in.
struct gssw_profile{
    // We keep one version, stored striped, for byte-sized alignment
	__m128i* profile_byte;	// 0: none
	// And another for word-sized alignment.
	__m128i* profile_word;	// 0: none
	const int8_t* read;
	const int8_t* mat;
	int32_t readLen;
	int32_t n;
	uint8_t bias;
};

//struct node;
//typedef struct node s_node;
typedef struct _gssw_node gssw_node;
typedef struct _gssw_node {
    void* data;
    uint32_t id;
    char* seq; // sequence
    int8_t* num; // numerical conversion of seq
    int32_t len; // length of sequence
    gssw_node** prev;
    int32_t count_prev;
    gssw_node** next;
    int32_t count_next;
    gssw_align* alignment;
} _gssw_node;

typedef struct {
    gssw_node* node;
    gssw_cigar cigar;
    int32_t ref_begin;
	int32_t ref_end;
	int32_t	read_begin;
	int32_t read_end;
} gssw_node_mapping;

typedef struct {
    gssw_node* node;
    gssw_alignment_end end;
} gssw_node_alignment_end;

typedef struct {
    uint32_t size;
    gssw_node* max_node;
    gssw_node** nodes;
} gssw_graph;

typedef struct {
    gssw_node* node;
    gssw_cigar* cigar;
} gssw_node_cigar;

typedef struct {
    uint32_t length;   // number of nodes traversed
    gssw_node_cigar* elements; // describes traceback
} gssw_graph_cigar;

typedef struct {
    int32_t position; // position in first node
    int16_t score;
    gssw_graph_cigar cigar;
} gssw_graph_mapping;

// matrices in the SSW algorithm
typedef enum gssw_matrix_t {Match, ReadGap, RefGap} gssw_matrix_t;

// position in the DP structure
typedef struct {
    int32_t read_pos;
    int32_t ref_pos;
    gssw_matrix_t from_matrix;
    gssw_matrix_t to_matrix;
    gssw_node* from_node;
    gssw_node* to_node;
} gssw_trace_back_deflection;

// an array of positions where the traceback follows a suboptimal alignment
typedef struct {
    int16_t score;
    gssw_trace_back_deflection* deflections; // first deflection should be the start coordinates
    int32_t num_deflections;
} gssw_alternate_alignment_ends;

// a linked list of top alternate alignments
typedef struct gssw_multi_align_stack_node {
    gssw_alternate_alignment_ends* alt_alignment;
    struct gssw_multi_align_stack_node* next;
    struct gssw_multi_align_stack_node* prev;
} gssw_multi_align_stack_node;

typedef struct {
    int32_t current_size;
    int32_t capacity;
    gssw_multi_align_stack_node* bottom_scoring;
    gssw_multi_align_stack_node* top_scoring;
} gssw_multi_align_stack;


#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/*!	@function	Create the query profile using the query sequence.
	@param	read	pointer to the query sequence; the query sequence needs to be numbers
	@param	readLen	length of the query sequence
	@param	mat	pointer to the substitution matrix; mat needs to be corresponding to the read sequence
	@param	n	the square root of the number of elements in mat (mat has n*n elements)
	@param	start_full_length_bonus	a bonus for aligning the full length of an alignment at the start of a query
	@param	end_full_length_bonus	a bonus for aligning the full length of an alignment at the end of a query
	@param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
						your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
	@return	pointer to the query profile structure
	@note	example for parameter read and mat:
			If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
			Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
			//A  C  G  T
			  2 -2 -2 -2 //A
			 -2  2 -2 -2 //C
			 -2 -2  2 -2 //G
			 -2 -2 -2  2 //T
			mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
*/
gssw_profile* gssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n,
                         int8_t start_full_length_bonus, int8_t end_full_length_bonus, const int8_t score_size);

/*!	@function	Create the quality-score adjusted query profile using the query sequence and its quality scores.
	@param	read        pointer to the query sequence; the query sequence needs to be numbers
	@param	readLen     length of the query sequence
	@param	adj_mat     pointer to the adjusted substitution matrix; mat needs to be corresponding to the read sequence
                            - see gssw_adjusted_qual_matrix and gssw_scaled_adjusted_qual_matrix
	@param	n           the square root of the number of elements in mat (mat has n*n elements)
	@param	start_full_length_bonus	a bonus for aligning the full length of an alignment at the start of a query
	@param	end_full_length_bonus	a bonus for aligning the full length of an alignment at the end of a query
	@param	score_size  estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 please set 0; if
                        your estimated best alignment score >= 255, please set 1; if you don't know, please set 2
	@return	pointer to the query profile structure
	@note	example for parameter read and mat:
			If the query sequence is: ACGTATC, the sequence that read points to can be: 1234142
			Then if the penalty for match is 2 and for mismatch is -2, the substitution matrix of parameter mat will be:
			//A  C  G  T
			  2 -2 -2 -2 //A
			 -2  2 -2 -2 //C
			 -2 -2  2 -2 //G
			 -2 -2 -2  2 //T
			mat is the pointer to the array {2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2, -2, -2, -2, -2, 2}
	@note   score_size option from gssw_init is implicitly set to 1 since adjusted scores should be scaled to as large
            of value as possible anyway to increase sensitivity
*/
gssw_profile* gssw_qual_adj_init (const int8_t* read, const int8_t* qual, const int32_t readLen, const int8_t* adj_mat,
                                  const int32_t n, int8_t start_full_length_bonus, int8_t end_full_length_bonus,
                                  const int8_t score_size);

/*!	@function	Release the memory allocated by function ssw_init.
	@param	p	pointer to the query profile structure
*/
void gssw_init_destroy (gssw_profile* p);

gssw_align* gssw_align_create(void);


// @function	ssw alignment.
/*!	@function	Do Striped Smith-Waterman alignment.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
				function ssw_init
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty
	@param	weight_gapE	the absolute value of gap extension penalty
	@param	flag	bitwise FLAG; (from high to low) bit 5: when setted as 1, function ssw_align will return the best alignment
					beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
					< filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
					cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
  					will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
 					setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
					the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
	@param	filters	score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
 					decription of the flag parameter for detailed usage.)
	@param	filterd	distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
					the decription of the flag parameter for detailed usage.)
	@param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
					readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
					return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
					alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
					largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
					picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
					reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
					score from the unmasked elements.
	@return	pointer to the alignment result structure
	@note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
			and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
			while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
			0-based coordinate.
*/
gssw_align* gssw_ssw_align (const gssw_profile* prof,
                            const int8_t* ref,
                            int32_t refLen,
                            const uint8_t weight_gapO,
                            const uint8_t weight_gapE,
                            const uint8_t flag,
                            const uint16_t filters,
                            const int32_t filterd,
                            const int32_t maskLen);

// @function	ssw fill.
/*!	@function	Do Striped Smith-Waterman alignment, but do not generate cigar.  Only return best score, filled matrix, and last E vector.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence needs to be numbers and corresponding to the mat parameter of
				function ssw_init
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty
	@param	weight_gapE	the absolute value of gap extension penalty
	@param	flag	bitwise FLAG; (from high to low) bit 5: when setted as 1, function ssw_align will return the best alignment
					beginning position; bit 6: when setted as 1, if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1
					< filterd), (whatever bit 5 is setted) the function will return the best alignment beginning position and
					cigar; bit 7: when setted as 1, if the best alignment score >= filters, (whatever bit 5 is setted) the function
  					will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever bit 5, 6 or 7 is
 					setted) the function will always return the best alignment beginning position and cigar. When flag == 0, only
					the optimal and sub-optimal scores and the optimal alignment ending position will be returned.
	@param	filters	score filter: when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filters will be used (Please check the
 					decription of the flag parameter for detailed usage.)
	@param	filterd	distance filter: when bit 6 of flag is setted as 1 and bit 8 is setted as 0, filterd will be used (Please check
					the decription of the flag parameter for detailed usage.)
	@param	maskLen	The distance between the optimal and suboptimal alignment ending position >= maskLen. We suggest to use
					readLen/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function will NOT
					return the suboptimal alignment information. Detailed description of maskLen: After locating the optimal
					alignment ending position, the suboptimal alignment score can be heuristically found by checking the second
					largest score in the array that contains the maximal score of each column of the SW matrix. In order to avoid
					picking the scores that belong to the alignments sharing the partial best alignment, SSW C library masks the
					reference loci nearby (mask length = maskLen) the best alignment ending position and locates the second largest
					score from the unmasked elements.
	@return	pointer to the alignment result structure
	@note	Whatever the parameter flag is setted, this function will at least return the optimal and sub-optimal alignment score,
			and the optimal alignment ending positions on target and query sequences. If both bit 6 and 7 of the flag are setted
			while bit 8 is not, the function will return cigar only when both criteria are fulfilled. All returned positions are
			0-based coordinate.
*/
gssw_align* gssw_ssw_fill (const gssw_profile* prof,
                           const int8_t* ref,
                           const int32_t refLen,
                           const uint8_t weight_gapO,
                           const uint8_t weight_gapE,
                           const int32_t maskLen,
                           gssw_seed* seed);


/*!	@function	Release the memory allocated by function ssw_align.
	@param	a	pointer to the alignment result structure
*/
void gssw_align_destroy (gssw_align* a);

/*!	@function	Release the memory allocated for matrices and pvE in s_align.
	@param	a	pointer to the alignment result structure
*/
void gssw_align_clear_matrix_and_seed (gssw_align* a);

/*! @function       Print score matrix, determines stride from result score
    @param refLen   Reference length.
    @param readLen  Read length.
    @param result   Alignment result.
*/
void gssw_print_score_matrix (const char* ref,
                              int32_t refLen,
                              const char* read,
                              int32_t readLen,
                              gssw_align* alignment,
                              FILE* out);

void gssw_graph_print(gssw_graph* graph);
void gssw_graph_print_stderr(gssw_graph* graph);

/*! @function         Trace back alignment across score matrix stored in alignment structure
    @param alignment  Alignment structure.
    @param end        Alignment ending position.
*/
    
    
gssw_cigar* gssw_alignment_trace_back (gssw_node* node,
                                       gssw_multi_align_stack* alt_alignment_stack,
                                       gssw_alternate_alignment_ends* alignment_deflections,
                                       int32_t* deflection_idx,
                                       int32_t final_traceback,
                                       uint16_t* score,
                                       int32_t* refEnd,
                                       int32_t* readEnd,
                                       int32_t* refGapFlag,
                                       int32_t* readGapFlag,
                                       const char* ref,
                                       int32_t refLen,
                                       const char* read,
                                       int8_t* qual_num,
                                       int32_t readLen,
                                       int8_t* nt_table,
                                       int8_t* score_matrix,
                                       uint8_t gap_open,
                                       uint8_t gap_extension,
                                       int8_t start_full_length_bonus,
                                       int8_t end_full_length_bonus);
    
gssw_cigar* gssw_alignment_trace_back_byte (gssw_node* node,
                                            gssw_multi_align_stack* alt_alignment_stack,
                                            gssw_alternate_alignment_ends* alignment_deflections,
                                            int32_t* deflection_idx,
                                            int32_t final_traceback,
                                            uint16_t* score,
                                            int32_t* refEnd,
                                            int32_t* readEnd,
                                            int32_t* refGapFlag,
                                            int32_t* readGapFlag,
                                            const char* ref,
                                            int32_t refLen,
                                            const char* read,
                                            int8_t* qual_num,
                                            int32_t readLen,
                                            int8_t* nt_table,
                                            int8_t* score_matrix,
                                            uint8_t gap_open,
                                            uint8_t gap_extension,
                                            int8_t start_full_length_bonus,
                                            int8_t end_full_length_bonus);

gssw_cigar* gssw_alignment_trace_back_word (gssw_node* node,
                                            gssw_multi_align_stack* alt_alignment_stack,
                                            gssw_alternate_alignment_ends* alignment_deflections,
                                            int32_t* deflection_idx,
                                            int32_t final_traceback,
                                            uint16_t* score,
                                            int32_t* refEnd,
                                            int32_t* readEnd,
                                            int32_t* refGapFlag,
                                            int32_t* readGapFlag,
                                            const char* ref,
                                            int32_t refLen,
                                            const char* read,
                                            int8_t* qual_num,
                                            int32_t readLen,
                                            int8_t* nt_table,
                                            int8_t* score_matrix,
                                            uint8_t gap_open,
                                            uint8_t gap_extension,
                                            int8_t start_full_length_bonus,
                                            int8_t end_full_length_bonus);

// Compute and return the traceback from a graph for which the alignment DP has been performed.
gssw_graph_mapping* gssw_graph_trace_back (gssw_graph* graph,
                                           const char* read,
                                           int32_t readLen,
                                           int8_t* nt_table,
                                           int8_t* score_matrix,
                                           uint8_t gap_open,
                                           uint8_t gap_extension,
                                           int8_t start_full_length_bonus,
                                           int8_t end_full_length_bonus);
    
gssw_graph_mapping* gssw_graph_trace_back_qual_adj (gssw_graph* graph,
                                                    const char* read,
                                                    const char* qual,
                                                    int32_t readLen,
                                                    int8_t* nt_table,
                                                    int8_t* adj_score_matrix,
                                                    uint8_t gap_open,
                                                    uint8_t gap_extension,
                                                    int8_t start_full_length_bonus,
                                                    int8_t end_full_length_bonus);

// Computes the traceback ending with the final character of the read aligned to the final character
// of a given node
gssw_graph_mapping* gssw_graph_trace_back_pinned (gssw_graph* graph,
                                                  gssw_node* pinned_node,
                                                  const char* read,
                                                  int32_t readLen,
                                                  int8_t* nt_table,
                                                  int8_t* score_matrix,
                                                  uint8_t gap_open,
                                                  uint8_t gap_extension,
                                                  int8_t start_full_length_bonus,
                                                  int8_t end_full_length_bonus);
    
gssw_graph_mapping* gssw_graph_trace_back_pinned_qual_adj (gssw_graph* graph,
                                                           gssw_node* pinned_node,
                                                           const char* read,
                                                           const char* qual,
                                                           int32_t readLen,
                                                           int8_t* nt_table,
                                                           int8_t* adj_score_matrix,
                                                           uint8_t gap_open,
                                                           uint8_t gap_extension,
                                                           int8_t start_full_length_bonus,
                                                           int8_t end_full_length_bonus);

// Computes an arbitrary number of highest scoring tracebacks ending with the final character of the
// read aligned to the final character of a given node
gssw_graph_mapping** gssw_graph_trace_back_pinned_multi (gssw_graph* graph,
                                                         gssw_node* pinned_node,
                                                         int32_t num_tracebacks,
                                                         const char* read,
                                                         int32_t readLen,
                                                         int8_t* nt_table,
                                                         int8_t* score_matrix,
                                                         uint8_t gap_open,
                                                         uint8_t gap_extension,
                                                         int8_t start_full_length_bonus,
                                                         int8_t end_full_length_bonus);

gssw_graph_mapping** gssw_graph_trace_back_pinned_qual_adj_multi (gssw_graph* graph,
                                                                  gssw_node* pinned_node,
                                                                  int32_t num_tracebacks,
                                                                  const char* read,
                                                                  const char* qual,
                                                                  int32_t readLen,
                                                                  int8_t* nt_table,
                                                                  int8_t* adj_score_matrix,
                                                                  uint8_t gap_open,
                                                                  uint8_t gap_extension,
                                                                  int8_t start_full_length_bonus,
                                                                  int8_t end_full_length_bonus);


    
/*! @function         Return 1 if the alignment is in 16/128bit (byte sized) or 0 if word-sized.
    @param alignment  Alignment structure.
*/
int gssw_is_byte (gssw_align* alignment);

/*! @function         Generate a traceback of the given alignment, using H, E, and F matrices
    @param alignment  Alignment structure.
    @param readPos    Starting position of alignment in reference.
    @param readPos    Starting position of alignment in read.
*/
//cigar* traceback (s_align* alignment, int32_t readPos, int32_t refPos);

void gssw_profile_destroy(gssw_profile* prof);
void gssw_seed_destroy(gssw_seed* seed);
gssw_seed* gssw_create_seed_byte(int32_t readLen, gssw_node** prev, int32_t count);
gssw_seed* gssw_create_seed_word(int32_t readLen, gssw_node** prev, int32_t count);

void gssw_cigar_push_back(gssw_cigar* c, char type, uint32_t length);
void gssw_cigar_push_front(gssw_cigar* c, char type, uint32_t length);
void gssw_reverse_cigar(gssw_cigar* c);
void gssw_print_cigar(gssw_cigar* c, FILE* out);
void gssw_cigar_destroy(gssw_cigar* c);

gssw_node* gssw_node_create(void* data,
                            const uint32_t id,
                            const char* seq,
                            const int8_t* nt_table,
                            const int8_t* score_matrix);
void gssw_node_destroy(gssw_node* n);
void gssw_node_add_prev(gssw_node* n, gssw_node* m);
void gssw_node_add_next(gssw_node* n, gssw_node* m);
void gssw_nodes_add_edge(gssw_node* n, gssw_node* m);
void gssw_node_del_prev(gssw_node* n, gssw_node* m);
void gssw_node_del_next(gssw_node* n, gssw_node* m);
void gssw_nodes_del_edge(gssw_node* n, gssw_node* m);
void gssw_node_replace_prev(gssw_node* n, gssw_node* m, gssw_node* p);
void gssw_node_replace_next(gssw_node* n, gssw_node* m, gssw_node* p);


gssw_node*
gssw_node_fill (gssw_node* node,
                const gssw_profile* prof,
                const uint8_t weight_gapO,
                const uint8_t weight_gapE,
                const int32_t maskLen,
                const gssw_seed* seed);

gssw_graph*
gssw_graph_fill (gssw_graph* graph,
                 const char* read_seq,
                 const int8_t* nt_table,
                 const int8_t* score_matrix,
                 const uint8_t weight_gapO,
                 const uint8_t weight_gapE,
                 const int8_t start_full_length_bonus,
                 const int8_t end_full_length_bonus,
                 const int32_t maskLen,
                 const int8_t score_size);

gssw_graph*
gssw_graph_fill_qual_adj(gssw_graph* graph,
                         const char* read_seq,
                         const char* read_qual,
                         const int8_t* nt_table,
                         const int8_t* adj_score_matrix,
                         const uint8_t weight_gapO,
                         const uint8_t weight_gapE,
                         const int8_t start_full_length_bonus,
                         const int8_t end_full_length_bonus,
                         const int32_t maskLen,
                         const int8_t score_size);
    
gssw_graph*
gssw_graph_fill_pinned (gssw_graph* graph,
                        const char* read_seq,
                        const int8_t* nt_table,
                        const int8_t* score_matrix,
                        const uint8_t weight_gapO,
                        const uint8_t weight_gapE,
                        const int8_t start_full_length_bonus,
                        const int8_t end_full_length_bonus,
                        const int32_t maskLen,
                        const int8_t score_size);

gssw_graph*
gssw_graph_fill_pinned_qual_adj(gssw_graph* graph,
                                const char* read_seq,
                                const char* read_qual,
                                const int8_t* nt_table,
                                const int8_t* adj_score_matrix,
                                const uint8_t weight_gapO,
                                const uint8_t weight_gapE,
                                const int8_t start_full_length_bonus,
                                const int8_t end_full_length_bonus,
                                const int32_t maskLen,
                                const int8_t score_size);
    
gssw_graph* gssw_graph_create(uint32_t size);
int32_t gssw_graph_add_node(gssw_graph* graph,
                            gssw_node* node);
void gssw_graph_clear(gssw_graph* graph);
void gssw_graph_destroy(gssw_graph* graph);
void gssw_graph_print_score_matrices(gssw_graph* graph,
                                     const char* read,
                                     int32_t readLen,
                                     FILE* out);

gssw_graph_mapping* gssw_graph_mapping_create(void);
void gssw_graph_mapping_destroy(gssw_graph_mapping* m);
gssw_graph_cigar* gssw_graph_cigar_create(void);
void gssw_graph_cigar_destroy(gssw_graph_cigar* g);
void gssw_print_graph_cigar(gssw_graph_cigar* g, FILE* out);
void gssw_print_graph_mapping(gssw_graph_mapping* gm, FILE* out);
//void gssw_graph_cigar_to_string(gssw_graph_cigar* g);
//void gssw_graph_mapping_to_string(gssw_graph_mapping* gm);

// some utility functions
int8_t* gssw_create_score_matrix(int32_t match, int32_t mismatch);
int8_t* gssw_create_nt_table(void);
int8_t* gssw_create_num(const char* seq,
                        const int32_t len,
                        const int8_t* nt_table);
int8_t* gssw_create_qual_num(const char* qual,
                             const int32_t len);
int8_t gssw_max_qual(const int8_t* qual,
                      const int32_t len);

    
/* Numerically computes the base of the logarithm in the log-odds interpretation of the scoring matrix */
double gssw_recover_log_base(const int8_t* score_matrix, const double* char_freqs, uint32_t alphabet_size, double tol);
/* Convenient wrapper for DNA matrices */
double gssw_dna_recover_log_base(int8_t match, int8_t mismatch, double gc_content, double tol);
    
// functions for adjusting alignments for base quality (see also gssw_qual_adj_init)

/*! @function   Create a base-quality adjusted scoring matrix to use across reads
 *  @param  max_qual         highest quality score to compute adjustments for
 *  @param  gap_open         gap open penalty
 *  @param  gap_extend       gap extend penalty
 *  @param  score_matrix     match score matrix
 *                              - should be of length alphabet_size * alphabet_size
 *                              - query sequences correspond to columns
 *  @param  char_freqs       frequency of characters in the alphabet
 *                              - should be of length alphabet_size
 *                              - should sum to 1
 *  @param  alphabet_size    number of characters in alphabet
 *  @param  tol              numerical tolerance for computing base of logarithm underlying log-odds scores
 *                              - recommended, 1e-14 to 1e-12
 *  @return                  pointer to the adjusted matrix
 *  @note   matrix is indexed by (qual_score) x (ref_char) x (query_char)
 */
int8_t* gssw_adjusted_qual_matrix(uint8_t max_qual, const int8_t* score_matrix, const double* char_freqs,
                                  uint32_t alphabet_size, double tol);

/*! @function   Create a scaled base-quality adjusted scoring matrix to use across reads (scaling improves sensitivity)
 *  @param  max_score        highest score to scale to
 *  @param  max_qual         highest quality score to compute adjustments for
 *  @param  gap_open_out     address of gap open penalty (will be modified)
 *  @param  gap_extend_out   address of gap extend penalty (will be modified)
 *  @param  score_matrix     match score matrix
 *                              - should be of length alphabet_size * alphabet_size
 *                              - query sequences correspond to columns
 *  @param  char_freqs       frequency of characters in the alphabet
 *                              - should be of length alphabet_size
 *                              - should sum to 1
 *  @param  alphabet_size    number of characters in alphabet
 *  @param  tol              numerical tolerance for computing base of logarithm underlying log-odds scores
 *                              - recommended, 1e-14 to 1e-12
 *  @return                  pointer to the adjusted matrix
 *  @note   The scores located at gap_open_out and gap_extend_out will be modified
 *  @note   matrix is indexed by (qual_score) x (ref_char) x (query_char)
 */
int8_t* gssw_scaled_adjusted_qual_matrix(int8_t max_score, uint8_t max_qual, int8_t* gap_open_out, int8_t* gap_extend_out,
                                         const int8_t* score_matrix, const double* char_freqs, uint32_t alphabet_size,
                                         double tol);

/* Creates a new set of matrices with 0 scores in the final row/column for ambiguous characters */
int8_t* gssw_add_ambiguous_char_to_adjusted_matrix(int8_t* adj_mat, uint8_t max_qual, uint32_t alphabet_size);
    
/* Wrapper for gssw_adjusted_qual_matrix using simple parameterization of char_freqs and score_matrix
   Adds row and column of 0s for N bases automatically */
int8_t* gssw_dna_scaled_adjusted_qual_matrix(int8_t max_score, uint8_t max_qual, int8_t* gap_open_out,
                                             int8_t* gap_extend_out, int8_t match_score, int8_t mismatch_score,
                                             double gc_content, double tol);

    
/* Initializes a stack to keep track of top alternate alignments (up to the capacity) in score sorted order */
gssw_multi_align_stack* gssw_new_multi_align_stack(int32_t capacity);
    
/* Destructor for the stack */
void gssw_delete_multi_align_stack(gssw_multi_align_stack* stack);

/* Constructor for a node on the stack */
gssw_multi_align_stack_node* gssw_new_multi_align_stack_node(gssw_alternate_alignment_ends* alignment_suffix, int16_t score,
                                                             int32_t read_pos, int32_t ref_pos, gssw_node* from_node,
                                                             gssw_node* to_node, gssw_matrix_t from_matrix, gssw_matrix_t to_matrix);

/* Destructor for a node on the stack */
void gssw_delete_multi_align_stack_node(gssw_multi_align_stack_node* stack_node);

/* Adds an alignment onto the stack if there is room or if it is higher scoring than any of the alignments currently in it */
void gssw_add_alignment(gssw_multi_align_stack* stack, gssw_alternate_alignment_ends* alignment_suffix, int16_t score,
                        int32_t read_pos, int32_t ref_pos, gssw_node* from_node, gssw_node* to_node, gssw_matrix_t from_matrix,
                        gssw_matrix_t to_matrix);
    
/* Retuns the lowest score of an alternate alignment in the stack */
int16_t gssw_min_alt_alignment_score(gssw_multi_align_stack* stack);
    
    
#ifdef __cplusplus
}
#endif	// __cplusplus

#endif	// SSW_H
