#include <iostream>
#include <vector>
#include <tuple>

class Matrix{
    std::vector<std::vector<double>> matrix;
public:
    Matrix();
    Matrix(std::size_t n, std::size_t m);
    explicit Matrix(std::vector<std::vector<double>> m) : matrix(std::move(m)) {}; //explicit + move - ???
    Matrix &operator=(std::vector<std::vector<double>> m) { matrix = std::move(m); return *this; }
    Matrix &operator=(const Matrix &rhs) = default; // = { matrix = rhs.matrix; return *this; }
    std::pair<std::size_t, std::size_t> get_size();
    std::vector<std::vector<double>> get_matrix();
    void print_matrix();
    double operator[](const std::pair<std::size_t, std::size_t>& Index);
    void transpose();
    Matrix transposed();
    std::vector<double> string(std::size_t i);
    std::vector<double> column(std::size_t i);
    Matrix submatrix(std::size_t x, std::size_t y);
    double determinant();
    Matrix adjugate();
    Matrix inverse();
};

Matrix::Matrix() {
    matrix.assign(0, std::vector<double> (0)); //???
}

Matrix::Matrix(std::size_t n, std::size_t m) {
    matrix.assign(n, std::vector<double> (m, 0));
}

std::pair<std::size_t, std::size_t> Matrix::get_size(){
    if(this->matrix.empty())
        return {0, 0};
    return {this->matrix.size(), this->matrix[0].size()};
}

std::vector<std::vector<double>> Matrix::get_matrix() {
    return this->matrix;
}

void Matrix::print_matrix() {
    std::pair<std::size_t, std::size_t> size = this->get_size();
    //std::cout << "Printing matrix with size: {" << size.first << ", " << size.second << "}: " << std::endl;
    for(std::size_t i = 0; i < size.first; i++){
        for(std::size_t j = 0; j < size.second; j++){
            std::cout << this->matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    //std::cout << "end of printing matrix" << std::endl;
}

double Matrix::operator[](const std::pair<std::size_t, std::size_t> &Index) {
    std::size_t i = Index.first;
    std::size_t j = Index.second;
    auto size = this->get_size();
    if((i >= 0 and i < size.first) and (j >= 0 and j < size.second)) {
        return matrix[Index.first][Index.second];
    }
    throw std::out_of_range("incorrect index");
}

Matrix Matrix::transposed() {
    std::vector<std::vector<double>> m = this->matrix;
    auto size = this->get_size();
    std::vector<std::vector<double>> tm (size.second, std::vector<double> (size.first));
    for(std::size_t i = 0; i < size.second; i++){
        for(std::size_t j = 0; j < size.first; j++){
            tm[i][j] = m[j][i];
        }
    }
    return Matrix(tm);
}

void Matrix::transpose() {
    this->matrix = this->transposed().matrix;
}

std::vector<double> Matrix::string(std::size_t i) {
    if (i >= 0 and i < this->get_size().first)
        return this->matrix[i];
    throw std::out_of_range("incorrect index");
    return std::vector<double> ();
}

std::vector<double> Matrix::column(std::size_t i) {
    auto size = this->get_size();
    if (i >= 0 and i < size.second) {
        std::vector<double> column(size.first);
        for (std::size_t j = 0; j < size.first; j++){
            column[j] = this->matrix[j][i];
        }
        return column;
    }
    throw std::out_of_range("incorrect index");
    return std::vector<double> ();
}

Matrix Matrix::submatrix(std::size_t x, std::size_t y) {
    auto size = this->get_size();
    if (x >= 0 and y >= 0 and x < size.first and y < size.second) {
        std::vector<std::vector<double>> _matrix = this->matrix;
        std::vector<std::vector<double>> submatrix;
        std::vector<double> tmp;
        for (std::size_t i = 0; i < size.first; i++){
            if(i != x){
                for(std::size_t j = 0; j < size.second; j++){
                    if(j != y)
                        tmp.push_back(_matrix[i][j]);
                }
                submatrix.push_back(tmp);
                tmp.clear();
            }
        }
        return Matrix(submatrix);
    }
    throw std::out_of_range("incorrect index");
}

double Matrix::determinant() {
    auto size = this->get_size();
    if (size.first == size.second){
        if (size.first == 2){
            return this->operator[]({0, 0}) * this->operator[]({1, 1})  \
                   - this->operator[]({1, 0}) * this->operator[]({0, 1});
        }
        if (size.first == 1){
            return this->operator[]({0, 0});
        }
        double sum = 0;
        std::size_t j = 1; // раскладываем по первой строке
        for (std::size_t i = 0; i < size.first; i++){
            sum += ((j + i) % 2 == 0 ? 1 : -1) * this->operator[]({j, i}) * this->submatrix(j, i).determinant();
        }
        return sum;
    }
    throw std::invalid_argument("determinant of a non square matrix");
}

Matrix Matrix::adjugate() {
    auto size = this->get_size();
    if (size.first == size.second){
        std::vector<std::vector<double>> adjugate_matrix(size.first, std::vector<double> (size.second));
        for (std::size_t i = 0; i < size.first; i++){
            for (std::size_t j = 0; j < size.second; j++){
                adjugate_matrix[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * this->submatrix(i, j).determinant();
            }
        }
        return Matrix(adjugate_matrix);
    }
    throw std::invalid_argument("adjugate of a non square matrix");
}

inline Matrix operator*(Matrix lhs, double rhs){
    std::vector<std::vector<double>> res(lhs.get_size().first, std::vector<double> (lhs.get_size().second));
    for (std::size_t i = 0; i < lhs.get_size().first; i++){
        for (std::size_t j = 0; j < lhs.get_size().second; j++){
            res[i][j] = lhs[{i, j}] * rhs;
        }
    }
    return Matrix(res);
}

Matrix Matrix::inverse() {
    auto size = this->get_size();
    if (size.first == size.second){
        if (this->determinant() != 0)
            return Matrix(this->adjugate().transposed() * (1/this->determinant()));
        throw std::invalid_argument("inverse of a zero determinant matrix");
    }
    throw std::invalid_argument("inverse of a non square matrix");
}

inline Matrix operator+(Matrix lhs, Matrix rhs){
    auto lhs_size = lhs.get_size();
    auto rhs_size = rhs.get_size();
    if(lhs_size == rhs_size){
         if(lhs_size.first + lhs_size.second == 0){
             return lhs;
         }
         std::vector<std::vector<double>> res_matrix = lhs.get_matrix(), rhs_matrix = rhs.get_matrix();
         for(std::size_t i = 0; i < lhs_size.first; i++){
             for(std::size_t j = 0; j < lhs_size.second; j++){
                 res_matrix[i][j] += rhs_matrix[i][j];
             }
         }
         return Matrix(res_matrix);
    }
    throw std::invalid_argument("sum of matrix with different sizes");
}

inline Matrix operator-(Matrix lhs, Matrix rhs){
    auto lhs_size = lhs.get_size();
    auto rhs_size = rhs.get_size();
    if(lhs_size == rhs_size){
        if(lhs_size.first + lhs_size.second == 0){
            return lhs;
        }
        std::vector<std::vector<double>> res_matrix = lhs.get_matrix(), rhs_matrix = rhs.get_matrix();
        for(std::size_t i = 0; i < lhs_size.first; i++){
            for(std::size_t j = 0; j < lhs_size.second; j++){
                res_matrix[i][j] -= rhs_matrix[i][j];
            }
        }
        return Matrix(res_matrix);
    }
    throw std::invalid_argument("sum of matrix with different sizes");
}

double cartesian_product(std::vector<double> lhs, std::vector<double> rhs){
    std::size_t lhs_size = lhs.size();
    std::size_t rhs_size = rhs.size();
    if(lhs_size == rhs_size){
        double sum = 0;
        for (std::size_t i = 0; i < lhs_size; i++){
            sum += lhs[i] * rhs[i];
        }
        return sum;
    }
    throw std::invalid_argument("decart_multiplication of vectors with different sizes");
}

inline Matrix operator*(Matrix lhs, Matrix rhs){
    auto lhs_size = lhs.get_size();
    auto rhs_size = rhs.get_size();
    if(lhs_size.second == rhs_size.first){
        std::vector<std::vector<double>> lhs_matrix = lhs.get_matrix(), rhs_matrix = rhs.get_matrix();
        std::vector<std::vector<double>> res_matrix (lhs_size.first, std::vector<double> (rhs_size.second));
        for(std::size_t i = 0; i < lhs_size.first; i++){
            for(std::size_t j = 0; j < rhs_size.second; j++){
                res_matrix[i][j] = cartesian_product(lhs.string(i), rhs.column(j));
            }
        }
        return Matrix(res_matrix);
    }
    throw std::invalid_argument("multiply of matrix with incorrect sizes");
}

int main() {
    int test_case;
    std::cout << "choose test case: 1 - test base functional, 2 - test inverse interface: " << std::endl;
    std::cin >> test_case;
    if(test_case == 1) {
        std::cout << "Created empty matrix with size: (2, 2): " << std::endl;
        Matrix m1(2, 2);
        m1.print_matrix();
        std::cout << "Created matrix with size: (2, 2) filled with value = 2: " << std::endl;
        std::vector<std::vector<double>> v1(2, std::vector<double>(2, 2));
        Matrix m2(v1);
        m2.print_matrix();
        std::cout << "Create matrix, which is sum of two previous: " << std::endl;
        Matrix m3 = m1 + m2;
        m3.print_matrix();
        std::pair<std::size_t, std::size_t> Index = {1, 0};
        std::cout << "Get element with definite index:" << std::endl;
        std::cout << "element[" << Index.first << "][" << Index.second << "] =  " << m3[Index] << std::endl;
        std::cout << "Create matrix, which is diff of the first two matrices: " << std::endl;
        Matrix m4 = m1 - m2;
        m4.print_matrix();
        std::pair<std::size_t, std::size_t> size = {2, 3};
        std::vector<std::vector<double>> matrix5(size.first, std::vector<double>(size.second));
        for (std::size_t i = 0; i < size.first; i++) {
            for (std::size_t j = 0; j < size.second; j++) {
                matrix5[i][j] = i*size.first + j;
            }
        }
        std::cout << "Create some new matrix: " << std::endl;
        Matrix m5(matrix5), m6(matrix5);
        m5.print_matrix();
        m6.transpose();
        std::cout << "Create transposed matrix: " << std::endl;
        m6.print_matrix();
        std::cout << "Multiply these two matrices: " << std::endl;
        Matrix m7 = m5 * m6;
        m7.print_matrix();
        int a = 3;
        std::cout << "Multiply the result matrix by the value = " << a << ':' << std::endl;
        Matrix m8 = m7 * a;
        m8.print_matrix();
        std::cout << "Multiply these two matrices in a different order: " << std::endl;
        Matrix m9 = m6 * m5;
        m9.print_matrix();
        std::cout << "Find the determinant of the previous matrix: " << std::endl;
        std::cout << "determinant = " << m9.determinant() << std::endl;
        std::cout << "Find the submatrix(1, 1) and it's determinant of the previous matrix: " << std::endl;
        m9.submatrix(1, 1).print_matrix();
        std::cout << "determinant = " << m9.submatrix(1, 1).determinant() << std::endl;
    }
    else if(test_case == 2){
        int n, m;
        std::cout << "Enter size of matrix n, m: " << std::endl;
        std::cin >> n >> m;
        std::cout << "Enter matrix " << n << 'x' << m << ": " << std::endl;
        std::vector<std::vector<double>> _matrix(n, std::vector<double>(m));
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < m; j++) {
                std::cin >> _matrix[i][j];
            }
        }
        Matrix matrix(_matrix);
        std::cout << "Matrix determinant: " << matrix.determinant() << std::endl;
        Matrix inverse_matrix = matrix.inverse();
        std::cout << "Inverse matrix: " << std::endl;
        inverse_matrix.print_matrix();
        std::cout << "Product of matrix and it's inverse: " << std::endl;
        Matrix product = matrix * inverse_matrix;
        product.print_matrix();
    }
    return 0;
}
