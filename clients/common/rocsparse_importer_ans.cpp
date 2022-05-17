#include "rocsparse_importer_ans.hpp"


rocsparse_importer_ans::rocsparse_importer_ans(const std::string& filename_) : m_filename(filename_)
{
    f = fopen(this->m_filename.c_str(), "r");
}
rocsparse_importer_ans::~rocsparse_importer_ans() { fclose(f); }


template<typename I>
rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_ell(I* M, I* N, I* ell_stride)
{
    if(!f)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell(params): cannot open file '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }

    // read number of rows
    std::size_t read_bytes = fread(M, sizeof(I), 1, f);
    read_bytes = fread(ell_stride, sizeof(I), 1, f);
    read_bytes = fread(N, sizeof(I), 1, f);
    //*N = *ell_cols;

    if (read_bytes == 0)
    { 
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell(params): read_bytes == 0 '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }

    m_ell_size = (*ell_stride) * (*N);

    if (!m_ell_size)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell: ell_size = 0 '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }
    return rocsparse_status_success;
}

template<typename T, typename I>
rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_ell(I* ell_ind, T* ell_val)
{
     if(!f)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell(vectors): cannot open file '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }

    std::size_t read_bytes = fread(ell_ind, sizeof(I), m_ell_size, f);
    read_bytes = fread(ell_val, sizeof(T), m_ell_size, f);

    if (read_bytes == 0)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell(vectors): read bytes == 0 '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }
    return rocsparse_status_success;
}

template<typename I>
rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_coo(I* coo_elements)
{
    if(!f)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_coo(elements): cannot open file '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }
    std::size_t read_bytes = fread(&m_coo_nelements, sizeof(I), 1, f);

    if (read_bytes == 0)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_coo(elements): read bytes == 0 '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }

    *coo_elements = m_coo_nelements;

    return rocsparse_status_success;
}

template<typename T, typename I>
rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_coo(I* coo_row_ind, I* coo_col_ind, T* coo_val)
{
    if(!f)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_ell(vectors): cannot open file '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }
    
    std::size_t read_bytes = fread(coo_row_ind, sizeof(I), m_coo_nelements, f);
    read_bytes = fread(coo_col_ind, sizeof(I), m_coo_nelements, f);
    read_bytes = fread(coo_val, sizeof(T), m_coo_nelements, f);

    if (read_bytes == 0)
    {
        std::cerr << "rocsparse_importer_ans::rocsparse_import_sparse_coo(vectors): read bytes == 0 '"
                  << this->m_filename << "' " << std::endl;
        return rocsparse_status_internal_error;
    }

    return rocsparse_status_success;
}


#define INSTANTIATE_TI(T, I)                                                                                                    \
    template rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_ell( I* ell_ind, T* ell_val);                     \
    template rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_coo( I* coo_row_ind, I* coo_col_ind, T* coo_val)
   
#define INSTANTIATE_I(I)                                                                                        \
    template rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_ell(I* M, I* N, I* ell_stride);   \
    template rocsparse_status rocsparse_importer_ans::rocsparse_import_sparse_coo(I* coo_elements)

INSTANTIATE_I(int32_t);
INSTANTIATE_I(int64_t);

INSTANTIATE_TI(double, int32_t);
INSTANTIATE_TI(double, int64_t);