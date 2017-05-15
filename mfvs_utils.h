/**
 * Copyright 2012 Alexandre Blondin Masse
 *
 * This file is part of the MFVS project.
 *
 * MFVS is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or any later version.

 * MFVS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>

using namespace std;

namespace mfvs {

/**
 * Breaks a string into tokens according to the
 * given separators.
 *
 * For instance, if the given string is
 * "jack:smith:28", the delimiter is ":", the
 * the tokens variable will contain the values
 * "jack", "smith" and "28" in a vector.
 *
 * @param str         the string to be broken
 *                    in tokens
 * @param tokens      the resulting tokens after
 *                    the break takes place
 * @param delimiters  the characters used as
 *                    separators
 */
void tokenize(const string& str, vector<string>& tokens,
              const string& delimiters = " ");

}