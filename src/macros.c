/* gb: Gröbner Basis
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gb.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file macros.c
 * \brief Macros for debugging and other stuff
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

/* note: idea of implementing debug output in this way comes from roman pearce
 *       and michael monagans maple mgb implementation */

/* debug macros */
#define SYMDBG  1 /* symbolic preprocessing */
#define SELDBG  1 /* spair selection */
#define MATDBG  1 /* matrix / linear algebra*/
#define UPDDBG  1 /* update process */
#define TIMDBG  1 /* timings */

#define DEBUG(s,...) do { if (s) printf(__VA_ARGS__); } while (0)
