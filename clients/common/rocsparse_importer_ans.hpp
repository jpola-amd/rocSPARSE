/* Importer ANSYS HYB FORMAT */

#pragma once
#ifndef ROCSPARSE_IMPORTER_ANS_HPP
#define ROCSPARSE_IMPORTER_ANS_HPP

#include "rocsparse_importer.hpp"


class rocsparse_importer_ans : public rocsparse_importer<rocsparse_importer_ans>
{
protected:
    std::string m_filename;
    int m_coo_nelements {};
    int m_ell_size {};

public:
    rocsparse_importer_ans(const std::string& filename_);
    ~rocsparse_importer_ans();

    //ASSUMPTION nnz = coo_val.size() + ell_val() (which is not entirely true)
    //M = nrows
    //N = ell_cols <---

    template<typename I = rocsparse_int>
    rocsparse_status rocsparse_import_sparse_ell(I* M, I* N, I* ell_stride);

    //WARNING T must be always double!!!
    template<typename T = double, typename I = rocsparse_int>
    rocsparse_status rocsparse_import_sparse_ell(I* ell_ind, T* ell_val);

    template<typename I>
    rocsparse_status rocsparse_import_sparse_coo(I* coo_elements);

    template<typename T = double, typename I = rocsparse_int>
    rocsparse_status rocsparse_import_sparse_coo(I* coo_row_ind, I* coo_col_ind, T* coo_val);

    


private:
    FILE* f;
    


};



#endif //ROCSPARSE_IMPORTER_ANS_HPP