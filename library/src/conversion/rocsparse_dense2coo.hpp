/*! \file */
/* ************************************************************************
* Copyright (c) 2020-2021 Advanced Micro Devices, Inc.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*
* ************************************************************************ */

#pragma once
#ifndef ROCSPARSE_DENSE2COO_HPP
#define ROCSPARSE_DENSE2COO_HPP

#include "handle.h"

template <typename I, typename T>
rocsparse_status rocsparse_dense2coo_template(rocsparse_handle          handle,
                                              rocsparse_order           order,
                                              I                         m,
                                              I                         n,
                                              const rocsparse_mat_descr descr,
                                              const T*                  A,
                                              I                         ld,
                                              const I*                  nnz_per_rows,
                                              T*                        coo_val,
                                              I*                        coo_row_ind,
                                              I*                        coo_col_ind);

#endif // ROCSPARSE_DENSE2COO_HPP
