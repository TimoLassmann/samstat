/*
 
 Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>
 
 This file is part of TagDust.
 
 TagDust is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 TagDust is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.
 
 */

/*! \file nuc_code.h
 \brief Declares arrays to convert nucleotide sequences. 
 
 Initializes nucleotide alphabet needed to parse input.
 \author Timo Lassmann
 \bug No known bugs.
 */

/**
 * @brief Converts nucleotides in ASCII to 0-5. 
 *
 *
 */
unsigned int nuc_code[256];

/**
 * @brief Converts 0-5 nucleotides into printable ASCII to.
 *
 *
 */
unsigned int rev_nuc_code[5];

void init_nuc_code(void);


