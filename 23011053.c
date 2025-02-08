#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
//AHMET FURKAN CINARKA, SAYISAL ANALIZ PRJ.
int menu(){
    int userchoice;
    printf("\nEnter key corresponding to the function:\n");
    printf("\n Negative Values are only Valued correctly When it's already in the string. Inputing negative x values may not work correctly. Trigonometric and Logarithmics are not supported.\n");
    printf("0-Exit\n1-Bisection\n2-Regula Falsi\n3-Newton Raphson\n4-Reverse Matrix\n5-Gauss Elimination\n");
    printf("6-Gauss Seidal\n7-Numerical Derivative\n8-Simpson Method\n9-Trapezoid Method\n10-Gregory Newton\n");
    printf("Your Choice:");
    scanf("%d",&userchoice);
    return userchoice;
}
// Will try to implement if i got the time
char *functionallocator() {
    // Allocate memory for 200 chars(Bigger functions may cause unknown instability issues)
    char *data = (char*)malloc(200 * sizeof(char));
    if (data == NULL) {
        // Memory allocation failure safeguard
        return NULL;
    }
    printf("\nf(x)=");
    // Consume newline character from the input buffer
    getchar();
    fgets((char*)data, 200, stdin);

    return (char*)data;
}
//String Parser comes Here ------
typedef struct {
    char *array;
    int top;
    int capacity;
} Stack;
// Function to create a new stack
Stack* createStack(int capacity) {
    Stack *stack = (Stack*)malloc(sizeof(Stack));
    stack->capacity = capacity;
    stack->top = -1;
    stack->array = (char*)malloc((size_t)stack->capacity * sizeof(char));
    return stack;
}
// Function to check if the stack is full
int isFull(Stack *stack) {
    return stack->top == stack->capacity - 1;
}
// Function to check if the stack is empty
int isEmpty(Stack *stack) {
    return stack->top == -1;
}
// Function to push an element onto the stack
void push(Stack *stack, char item) {
    if (isFull(stack)) {
        printf("Stack Overflow\n");
        return;
    }
    stack->array[++stack->top] = item;
}
// Function to pop an element from the stack
char pop(Stack *stack) {
    if (isEmpty(stack)) {
        printf("Stack Underflow\n");
        return '\0';
    }
    return stack->array[stack->top--];
}
// Function to get the top element of the stack without popping it
char peek(Stack *stack) {
    if (isEmpty(stack)) {
        printf("Stack is empty\n");
        return '\0';
    }
    return stack->array[stack->top];
}
// Function to determine the precedence of operators
int precedence(char op) {
    if (op == '^')
        return 3;
    else if (op == '*' || op == '/')
        return 2;
    else if (op == '+' || op == '-')
        return 1;
    return 0;
}
char* modifyString(char *strupdate, double xvalue) {
    int len = (int)strlen(strupdate);
    int index = 0;

    // Check if the first character is '-'
    if (strupdate[index] == '-') {
        // If so, append a '0' before the '-' to make it '-0'
        memmove(&strupdate[index + 1], &strupdate[index], (size_t)(len - index + 1));
        strupdate[index] = '0';
        index++;
        len++; // Increment length to account for the added '0'
    }

    while (index < len) {
        if (isdigit(strupdate[index]) && strupdate[index + 1] == 'x') {
            // If a number followed by 'x' is found
            // Replace 'x' with the value provided by the user
            strupdate[index + 1] = '*';
            char replacement[20]; // Temporary buffer for replacement number
            sprintf(replacement, "%.3lf", xvalue); // Convert float to string
            int replacementLen = (int)strlen(replacement);
            // Shift the rest of the string to accommodate the replacement number
            memmove(&strupdate[index + 2 + replacementLen], &strupdate[index + 2], (size_t)(len - index - 1));
            // Copy the replacement number into the string
            memcpy(&strupdate[index + 2], replacement, (size_t)replacementLen);
            len += replacementLen - 1; // Update the length of the string
            index += replacementLen + 1; // Move to the next non-modified character
        } else if (strupdate[index] == 'x' && (index == 0 || !isdigit(strupdate[index - 1]))) {
            // If 'x' is found without a preceding number, replace with given x value
            char replacement[20]; // Temporary buffer for replacement number
            sprintf(replacement, "%.3lf", xvalue); // Convert float to string
            int replacementLen = (int)strlen(replacement);
            // Shift the rest of the string to accommodate the replacement number
            memmove(&strupdate[index + replacementLen], &strupdate[index + 1], (size_t)(len - index));
            // Copy the replacement number into the string
            memcpy(&strupdate[index], replacement, (size_t)replacementLen);
            len += replacementLen - 1; // Update the length of the string
            index += replacementLen; // Move to the next non-modified character
            
        } else if (strupdate[index]=='e'){
        	strupdate[index]='3';
		}
		
		
		 else {
            index++; // Move to the next character
        }
    }

    strupdate[len+2] = '\0'; // Null-terminate the modified string
    return strupdate;
}
//INFIX ------
// Function to convert infix expression to postfix
char* infixToPostfix(char* expression) {
    int len = (int)strlen(expression);
    Stack *stack = createStack(len); // Stack to store operators
    char *postfix = (char*)malloc((size_t)(2 * len + 1) * sizeof(char)); // Output postfix expression
    int k = 0; // Index for output postfix expression

    int i;for (i = 0; i < len; i++) {
        char ch = expression[i];
        if (isdigit(ch) || ch == '.') {
            int j = i;
            while (isdigit(expression[j]) || expression[j] == '.') {
                postfix[k++] = expression[j++];
            }
            postfix[k++] = ' '; // Add space after operand
            i = j - 1; // Update i to the last digit of the operand
        } else if (isalnum(ch)) {
            postfix[k++] = ch;
            postfix[k++] = ' '; // Add space after operand
        } else if (ch == '(') {
            push(stack, ch);
        } else if (ch == ')') {
            while (!isEmpty(stack) && peek(stack) != '(') {
                postfix[k++] = pop(stack);
                postfix[k++] = ' '; // Add space after operand
            }
            if (!isEmpty(stack) && peek(stack) != '(') {
                printf("Invalid expression\n");
                return NULL;
            } else {
                pop(stack); // Discard the '('
            }
        } else { // Operator
            while (!isEmpty(stack) && precedence(ch) <= precedence(peek(stack))) {
                postfix[k++] = pop(stack);
                postfix[k++] = ' '; // Add space after operand
            }
            push(stack, ch);
        }
    }

    // Pop remaining operators from the stack
    while (!isEmpty(stack)) {
        postfix[k++] = pop(stack);
        postfix[k++] = ' '; // Add space after operand
    }

    postfix[k] = '\0'; // Null-terminate the postfix expression
    return postfix;
}
//POSTFIX ------
// Stack structure for operands
typedef struct {
    float *array;
    int top;
    int capacity;
} OperandStack;
// Function to create a new operand stack
OperandStack* createOperandStack(int capacity) {
    OperandStack *stack = (OperandStack*)malloc(sizeof(OperandStack));
    stack->capacity = capacity;
    stack->top = -1;
    stack->array = (float*)malloc((size_t)stack->capacity * sizeof(float));
    return stack;
}
// Function to check if the operand stack is full
int isOperandStackFull(OperandStack *stack) {
    return stack->top == stack->capacity - 1;
}
// Function to check if the operand stack is empty
int isOperandStackEmpty(OperandStack *stack) {
    return stack->top == -1;
}
// Function to push an operand onto the stack
void pushOperand(OperandStack *stack, float operand) {
    if (isOperandStackFull(stack)) {
        printf("Operand Stack Overflow\n");
        return;
    }
    stack->array[++stack->top] = operand;
}
// Function to pop an operand from the stack
float popOperand(OperandStack *stack) {
    if (isOperandStackEmpty(stack)) {
        printf("Operand Stack Underflow\n");
        return -1; // Error value
    }
    return stack->array[stack->top--];
}
float parseOperand(char* expression, int* index) {
    float operand = 0;
    int decimalFound = 0;
    float decimalPlace = 0.1f; // For handling decimal places

    while (isdigit(expression[*index]) || expression[*index] == '.') {
        if (expression[*index] == '.') {
            decimalFound = 1;
            (*index)++;
            continue;
        }
        if (!decimalFound) {
            operand = operand * 10 + (float)(expression[*index] - '0');
        } else {
            operand = operand + (float)(expression[*index] - '0') * decimalPlace;
            decimalPlace /= 10;
        }
        (*index)++;
    }
    return operand;
}
float evaluatePostfix(char* expression) {
    OperandStack *stack = createOperandStack(100); // Assuming maximum expression length is 100
    int len = (int)strlen(expression);

    int i;for (i = 0; i < len; i++) {
        char ch = expression[i];
        printf("Processing character: %c\n", ch);
        if (isdigit(ch) || ch == '.') {
            float operand = parseOperand(expression, &i);
            printf("Pushing operand: %.3f\n", operand);
            pushOperand(stack, operand);
        } else if (ch == '+' || ch == '-' || ch == '*' || ch == '/') {
            if (isOperandStackEmpty(stack)) {
                printf("Operand Stack Underflow\n");
                free(stack->array);
                free(stack);
                return -1; // Error value
            }
            float operand2 = popOperand(stack);
            if (isOperandStackEmpty(stack)) {
                printf("Operand Stack Underflow\n");
                free(stack->array);
                free(stack);
                return -1; // Error value
            }
            float operand1 = popOperand(stack);
            float result;
            printf("Popped operands: %.3f, %.3f\n", operand1, operand2);
            switch (ch) {
                case '+':
                    result = operand1 + operand2;
                    break;
                case '-':
                    result = operand1 - operand2;
                    break;
                case '*':
                    result = operand1 * operand2;
                    break;
                case '/':
                    if (operand2 == 0) {
                        printf("Division by zero error\n");
                        free(stack->array);
                        free(stack);
                        return -1; // Error value
                    }
                    result = operand1 / operand2;
                    break;
                default:
                    printf("Invalid operator\n");
                    free(stack->array);
                    free(stack);
                    return -1; // Error value
            }
            printf("Result of operation: %.3f\n", result);
            pushOperand(stack, result);
        } else if (ch == '^') {
            if (isOperandStackEmpty(stack)) {
                printf("Operand Stack Underflow\n");
                free(stack->array);
                free(stack);
                return -1; // Error value
            }
            float exponent = popOperand(stack);
            if (isOperandStackEmpty(stack)) {
                printf("Operand Stack Underflow\n");
                free(stack->array);
                free(stack);
                return -1; // Error value
            }
            float base = popOperand(stack);
            float result = (float)pow(base, exponent);
            printf("Result of operation: %.3f\n", result);
            pushOperand(stack, result);
        }
    }

    if (isOperandStackEmpty(stack)) {
        printf("Operand Stack Underflow\n");
        free(stack->array);
        free(stack);
        return -1; // Error value
    } else {
        float finalResult = popOperand(stack);
        printf("Final result: %.3f\n", finalResult);
        free(stack->array);
        free(stack);
        return finalResult;
    }
}
//String Parser ends here------
double stringparsercall(double x, char* datap) {
    double xvalue;
    double result;
    xvalue = x;

    int size = (int)strlen(datap) + 2;

    char *data;
    data = (char*)malloc((size_t)size);

    memcpy(data, datap, (size_t)size);
    
    char *modifiedStr = modifyString(data, xvalue);
    printf("Infix expression: %s", modifiedStr);
    char *postfixExpression = infixToPostfix(modifiedStr);
    printf("\nPostfix expression: %s\n", postfixExpression);
    result = evaluatePostfix(postfixExpression);
    printf("\n %.3lf is the result.\n", result);

    free(modifiedStr);

    return result;
}
// Bisection method function
double bisection(double a, double b, double tol, char * functiondata) {
    double c, fa, fb, fc;

    fa = stringparsercall(a, functiondata);
    fb = stringparsercall(b, functiondata);
    if (fa * fb >= 0) {
        printf("Function does not guarantee a root within [%lf, %lf]\n", a, b);
        return NAN;
    }

    while ((b - a) / 2.0 > tol) {
        c = (a + b) / 2.0;
        fc = stringparsercall(c, functiondata);
        
        printf("a=%lf, fa=%lf, b=%lf, fb=%lf, c=%lf, fc=%lf\n", a, fa, b, fb, c, fc);
        
        if (fc == 0.0 || (b - a) / 2.0 < tol) {
            printf("%lf is the root.\n", c);
            return c;
        }
        
        if (fc * fa < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    c = (a + b) / 2.0;  // Final midpoint is our best estimate for root
    printf("%lf is the approximate root.\n", c);
    return c;
}
double regulaFalsi(double a, double b, char* functiondata, double tol, int max_iter) {
    double fa = stringparsercall(a, functiondata);
    double fb = stringparsercall(b, functiondata);

    if (fa * fb >= 0) {
        printf("The function must have different signs at a and b.\n");
        return NAN;  // Error indicator
    }

    double c = a;  // Initialize c
    int i;for (i = 0; i < max_iter; i++) {
        // Calculate the point where the function changes sign
        c = (a * fb - b * fa) / (fb - fa);
        double fc = stringparsercall(c, functiondata);

        // Check if we found the root or if the tolerance is met
        if (fabs(fc) < tol) {
            //used fabs to get absolute value of fc. 
            return c;
        }

        // Decide the side to repeat the steps
        if (fc * fa < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        
        printf("Iteration %d: a=%lf, b=%lf, c=%lf, fa=%lf, fb=%lf, fc=%lf\n", i+1, a, b, c, fa, fb, fc);
    }

    printf("The method did not converge after %d iterations.\n", max_iter);
    return c;  // Return the best estimate of the root
}
double newtonRaphson(double initial_guess, char* functiondata, double tol, int max_iter) {
    double x = initial_guess;  // Starting value for root finding
    double h = 0.001;  // Smaller step size for derivative approximation
    int iter = 0;
    double fx, dfx, x_next;

    do {
        fx = stringparsercall(x, functiondata);  // Calculate f(x)
        // Calculate derivative f'(x) using central difference formula
        dfx = (stringparsercall(x + h, functiondata) - stringparsercall(x - h, functiondata)) / (2 * h);

        if (fabs(dfx) < tol) {  // Check if derivative is too small
            printf("Derivative too small.\n");
            return NAN;
        }

        x_next = x - fx / dfx;  // Newton-Raphson formula

        printf("Iteration %d: x = %lf, f(x) = %lf, f'(x) = %lf, x_next = %lf\n", iter, x, fx, dfx, x_next);

        if (fabs(x_next - x) < tol) {  // Check convergence
            printf("Root found: %lf\n", x_next);
            return x_next;
        }

        x = x_next;  // Update x to the new value
        iter++;
    } while (iter < max_iter);

    printf("Method did not converge after %d iterations.\n", max_iter);
    return x_next;  // Return the last computed value of x
}
//Inverse Matrix Function prototypes
float calculateDeterminant(float **matrix, int size);
void cofactorMatrix(float **matrix, int size);
void transposeMatrix(float **matrix, float **result, int size);
void freeMatrix(float **matrix, int size);
//Inverse Matrix functions
float calculateDeterminant(float **matrix, int size) {
    float sign = 1, det = 0, **submatrix;
    int i, j, m, n, c;

    if (size == 1)
        return (matrix[0][0]);
    else {
        det = 0;
        submatrix = (float **)malloc((size_t)(size - 1) * sizeof(float *));
        if (submatrix == NULL) {
            printf("Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < size - 1; i++) {
            submatrix[i] = (float *)malloc((size_t)(size - 1) * sizeof(float));
            if (submatrix[i] == NULL) {
                printf("Memory allocation failed\n");
                int k;for (k = 0; k < i; k++) {
                    free(submatrix[k]);
                }
                free(submatrix);
                exit(EXIT_FAILURE);
            }
        }

        for (c = 0; c < size; c++) {
            m = 0;
            n = 0;

            for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    if (i != 0 && j != c) {
                        submatrix[m][n] = matrix[i][j];

                        if (n < (size - 2))
                            n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }

            det = det + sign * (matrix[0][c] * calculateDeterminant(submatrix, size - 1));
            sign = -1 * sign;
        }

        for (i = 0; i < size - 1; i++) {
            free(submatrix[i]);
        }
        free(submatrix);
    }

    return (det);
}
void cofactorMatrix(float **matrix, int size) {
    float **minorMatrix, **cofactor;
    int row, col, subrow, subcol, i, j;

    minorMatrix = (float **)malloc((size_t)(size - 1) * sizeof(float *));
    if (minorMatrix == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size - 1; i++) {
        minorMatrix[i] = (float *)malloc((size_t)(size - 1) * sizeof(float));
        if (minorMatrix[i] == NULL) {
            printf("Memory allocation failed\n");
            int k; for (k = 0; k < i; k++) {
                free(minorMatrix[k]);
            }
            free(minorMatrix);
            exit(EXIT_FAILURE);
        }
    }

    cofactor = (float **)malloc((size_t)size * sizeof(float *));
    if (cofactor == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; i++) {
        cofactor[i] = (float *)malloc((size_t)size * sizeof(float));
        if (cofactor[i] == NULL) {
            printf("Memory allocation failed\n");
            int k;for (k = 0; k < i; k++) {
                free(cofactor[k]);
            }
            free(cofactor);
            exit(EXIT_FAILURE);
        }
    }

    for (row = 0; row < size; row++) {
        for (col = 0; col < size; col++) {
            subrow = 0;
            subcol = 0;

            for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    if (i != row && j != col) {
                        minorMatrix[subrow][subcol] = matrix[i][j];

                        if (subcol < (size - 2))
                            subcol++;
                        else {
                            subcol = 0;
                            subrow++;
                        }
                    }
                }
            }

            cofactor[row][col] = (float)pow(-1, row + col) * calculateDeterminant(minorMatrix, size - 1);
        }
    }

    transposeMatrix(matrix, cofactor, size);

    for (i = 0; i < size - 1; i++) {
        free(minorMatrix[i]);
    }
    free(minorMatrix);

    for (i = 0; i < size; i++) {
        free(cofactor[i]);
    }
    free(cofactor);
}
void transposeMatrix(float **matrix, float **result, int size) {
    float **transpose, **inverse;
    int i, j;

    transpose = (float **)malloc((size_t)size * sizeof(float *));
    if (transpose == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; i++) {
        transpose[i] = (float *)malloc((size_t)size * sizeof(float));
        if (transpose[i] == NULL) {
            printf("Memory allocation failed\n");
            int k;for (k = 0; k < i; k++) {
                free(transpose[k]);
            }
            free(transpose);
            exit(EXIT_FAILURE);
        }
    }

    inverse = (float **)malloc((size_t)size * sizeof(float *));
    if (inverse == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < size; i++) {
        inverse[i] = (float *)malloc((size_t)size * sizeof(float));
        if (inverse[i] == NULL) {
            printf("Memory allocation failed\n");
            int k;for (k = 0; k < i; k++) {
                free(inverse[k]);
            }
            free(inverse);
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            transpose[i][j] = result[j][i];
        }
    }

    float det = calculateDeterminant(matrix, size);

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            inverse[i][j] = transpose[i][j] / det;
        }
    }

    printf("\n\n\nThe inverse of matrix is: \n");

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            printf("\t%f", inverse[i][j]);
        }
        printf("\n");
    }

    for (i = 0; i < size; i++) {
        free(transpose[i]);
        free(inverse[i]);
    }
    free(transpose);
    free(inverse);
}
void freeMatrix(float **matrix, int size) {
    int i;for (i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}
void inversematrixcall(){
    float **inputMatrix;
    float determinant;
    int order;

    printf("Enter the order of the Matrix: ");
    scanf("%d", &order);

    inputMatrix = (float **)malloc((size_t)order * sizeof(float *));
    if (inputMatrix == NULL) {
        printf("Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    int i;for (i = 0; i < order; i++) {
        inputMatrix[i] = (float *)malloc((size_t)order * sizeof(float));
        if (inputMatrix[i] == NULL) {
            printf("Memory allocation failed\n");
            int k;for (k = 0; k < i; k++) {
                free(inputMatrix[k]);
            }
            free(inputMatrix);
            exit(EXIT_FAILURE);
        }
    }

    printf("Enter the elements of %dX%d Matrix: \n", order, order);
    for (i = 0; i < order; i++) {
        int j;for (j = 0; j < order; j++) {
            printf("Enter [%d][%d]",i,j);
            scanf("%f", &inputMatrix[i][j]);
        }
    }

    determinant = calculateDeterminant(inputMatrix, order);

    if (determinant == 0)
        printf("\nInverse of Entered Matrix is not possible\n");
    else {
        cofactorMatrix(inputMatrix, order);
    }

    freeMatrix(inputMatrix, order);
}
//Gauss Elimination Method
void gausseliminationcall(){
    int row, col, k, order; //iterators for input and order flag
    float **matrix, *solutions, constant, sum=0.0f;

    printf("\nOrder of Matrix:");
    scanf("%d", &order);

    // Memory allocation for matrix:
    matrix = (float **)malloc((size_t)(order + 1) * sizeof(float *));
    for(row = 1; row <= order; row++) {
        matrix[row] = (float *)malloc((size_t)(order + 2) * sizeof(float));
    }

    // Memory allocation for solutions:
    solutions = (float *)malloc((size_t)(order + 1) * sizeof(float));

    printf("\nEnter the elements of the augmented matrix row-wise:\n\n");
    for(row = 1; row <= order; row++) {
        for(col = 1; col <= (order + 1); col++) {
            printf("Matrix[%d][%d] : ", row, col);
            scanf("%f", &matrix[row][col]);
        }
    }

    for(col = 1; col <= order; col++) {
        for(row = 1; row <= order; row++) {
            if(row > col) {
                constant = matrix[row][col] / matrix[col][col];
                for(k = 1; k <= order + 1; k++) {
                    matrix[row][k] = matrix[row][k] - constant * matrix[col][k];
                }
            }
        }
    }

    solutions[order] = matrix[order][order + 1] / matrix[order][order];

    for(row = order - 1; row >= 1; row--) {
        sum = 0;
        for(col = row + 1; col <= order; col++) {
            sum = sum + matrix[row][col] * solutions[col];
        }
        solutions[row] = (matrix[row][order + 1] - sum) / matrix[row][row];
    }

    printf("\nThe solution is: \n");
    for(row = 1; row <= order; row++) {
        printf("\nx%d=%f\t", row, solutions[row]);
    }

    // Freeing memory to ensure no memory leaks etc:
    for(row = 1; row <= order; row++) {
        free(matrix[row]);
    }
    free(matrix);
    free(solutions);
}
void gaussSeidel(float **a, float *b, float *x, int n, float tol, int max_iter) {
	int i,j;
    float *x_old = malloc(n * sizeof(float));
    if (x_old == NULL) {
        printf("Memory allocation failed.\n");
        return;
    }

    for (i = 0; i < n; i++) {
        x[i] = 0.0; // Initial guess
    }

    int iter = 0;
    float norm;
    do {
        for (i = 0; i < n; i++) {
            x_old[i] = x[i];
        }

        for (i = 0; i < n; i++) {
            float sum = b[i];
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sum -= a[i][j] * x[j];
                }
            }
            x[i] = sum / a[i][i];

            if (a[i][i] == 0) {
                printf("Division by zero.\n");
                free(x_old);
                return;
            }
        }

        norm = 0.0;
        for (i = 0; i < n; i++) {
            norm += fabs(x[i] - x_old[i]);
        }

        iter++;
    } while (norm > tol && iter < max_iter);

    if (iter == max_iter) {
        printf("Gauss-Seidel method did not converge after %d iterations.\n", max_iter);
    } else {
        printf("Converged after %d iterations.\n", iter);
    }

    free(x_old);
}
void gaussSeidelCall() {
    int n,i,k,j;
    printf("Enter the number of variables: ");
    scanf("%d", &n);

    float **a = malloc(n * sizeof(float*));
    float *b = malloc(n * sizeof(float));
    float *x = malloc(n * sizeof(float));
    
    if (a == NULL || b == NULL || x == NULL) {
        printf("Memory allocation failed.\n");
        return;
    }

    for (i = 0; i < n; i++) {
        a[i] = malloc(n * sizeof(float));
        if (a[i] == NULL) {
            printf("Memory allocation failed.\n");
            for (k = 0; k < i; k++) {
                free(a[k]);
            }
            free(a);
            free(b);
            free(x);
            return;
        }
    }

    printf("Enter the coefficients of the matrix (row by row):\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            scanf("%f", &a[i][j]);
        }
    }

    printf("Enter the constants vector (b):\n");
    for (i = 0; i < n; i++) {
        scanf("%f", &b[i]);
    }

    float tol;
    int max_iter;
    printf("Enter tolerance and maximum iterations:\n");
    scanf("%f %d", &tol, &max_iter);

    gaussSeidel(a, b, x, n, tol, max_iter);

    printf("The solution is:\n");
    for (i = 0; i < n; i++) {
        printf("x%d = %f\n", i + 1, x[i]);
    }

    for (i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
}
// Numerical Derivative Calculation
double numericalDerivative(double x, char *functiondata, double h) {
    double fxh1 = stringparsercall(x + h, functiondata);  // f(x + h)
    double fxh2 = stringparsercall(x - h, functiondata);  // f(x - h)
    double derivative = (fxh1 - fxh2) / (2 * h);          // Central difference
    return derivative;
}
void numericalDerivativeCall() {
    char *functiondata = functionallocator();
    double x, h;
    printf("Enter the point x at which derivative is calculated: ");
    scanf("%lf", &x);
    printf("Enter the step size h: ");
    scanf("%lf", &h);
    double derivative = numericalDerivative(x, functiondata, h);
    printf("The derivative of the function at x = %lf is %lf\n", x, derivative);
    free(functiondata);
}
// Simpson's 1/3 Method
double simpsonOneThird(double a, double b, int n, char *functiondata) {
	int i;
    if (n % 2 != 0) {
        printf("Error: n must be even.\n");
        return NAN; // Simpson's 1/3 rule requires an even number of intervals.
    }

    double h = (b - a) / n;
    double sum = stringparsercall(a, functiondata) + stringparsercall(b, functiondata); // f(x_0) and f(x_n)
    
    for (i = 1; i < n; i++) {
        double x = a + i * h;
        double fx = stringparsercall(x, functiondata);
        if (i % 2 == 0) { // Even index terms have a coefficient of 2
            sum += 2 * fx;
        } else { // Odd index terms have a coefficient of 4
            sum += 4 * fx;
        }
    }

    double result = (h / 3) * sum;
    return result;
}
void simpsonOneThirdCall() {
    char *functiondata = functionallocator();
    double a, b;
    int n;

    printf("Enter the lower limit a: ");
    scanf("%lf", &a);
    printf("Enter the upper limit b: ");
    scanf("%lf", &b);
    printf("Enter the number of subintervals n (must be even): ");
    scanf("%d", &n);

    double result = simpsonOneThird(a, b, n, functiondata);
    printf("The approximate value of the integral is %lf\n", result);

    free(functiondata);
}
// Trapezoidal Rule for numerical integration
double trapezoidalRule(double a, double b, int n, char *functiondata) {
	int i;
    if (n < 1) {
        printf("Error: n must be at least 1.\n");
        return NAN;
    }

    double h = (b - a) / n;
    double sum = stringparsercall(a, functiondata);  // f(a)
    double x;

    for (i = 1; i < n; i++) {
        x = a + i * h;
        sum += 2 * stringparsercall(x, functiondata);  // 2 * sum of f(x_i)
    }

    sum += stringparsercall(b, functiondata);  // f(b)
    return (h / 2) * sum;
}
void trapezoidalRuleCall() {
    char *functiondata = functionallocator();
    double a, b;
    int n;

    printf("Enter the lower limit a: ");
    scanf("%lf", &a);
    printf("Enter the upper limit b: ");
    scanf("%lf", &b);
    printf("Enter the number of subintervals n: ");
    scanf("%d", &n);

    double result = trapezoidalRule(a, b, n, functiondata);
    printf("The approximate value of the integral using the Trapezoidal Rule is %lf\n", result);

    free(functiondata);
}
// Function to calculate factorial (i need it in gregory newton)
int factorial(int n) {
    if (n == 0 || n == 1) return 1;
    else return n * factorial(n - 1);
}
// Function to calculate the k-th backward difference
double backwardDifference(double *y, int n, int k) {
	int i,j;
    double diff[n];
    for (i = 0; i < n; i++) {
        diff[i] = y[i];
    }
    
    for (j = 1; j <= k; j++) {
        for (i = n-1; i >= j; i--) {
            diff[i] = diff[i] - diff[i-1];
        }
    }
    
    return diff[n-1]; // Return the last element of the difference for k-th backward
}
// Gregory-Newton Interpolation
double gregoryNewtonBackward(double *x, double *y, int n, double value) {
	int k;
    double h = x[1] - x[0]; // Calculate interval size
    double s = (value - x[n-1]) / h;
    double result = y[n-1]; // Start with the last y value
    double s_term = 1;
    
    for (k = 1; k < n; k++) {
        s_term *= (s + k - 1) / k;
        result += s_term * backwardDifference(y, n, k) / factorial(k);
    }

    return result;
}
void gregoryNewtonBackwardCall() {
    int n,i;
    printf("Enter the number of data points: ");
    scanf("%d", &n);

    double *x = malloc(n * sizeof(double));
    double *y = malloc(n * sizeof(double));

    printf("Enter the data points (x and y pairs):\n");
    for (i = 0; i < n; i++) {
        scanf("%lf %lf", &x[i], &y[i]);
    }

    double value;
    printf("Enter the x value for which you want to interpolate: ");
    scanf("%lf", &value);

    double result = gregoryNewtonBackward(x, y, n, value);
    printf("The interpolated value at x = %lf is %lf\n", value, result);

    free(x);
    free(y);
}
int methodchooser(int userchoice){
    char *functiondata;
    double equationresult;
    double a, b, c;
    int max;
    switch(userchoice){
        case 1:
        	//DONE
            // Bisection
            printf("\nBisection:\n");
            functiondata = functionallocator();
            printf("Enter lowerbound, upperbound, tolerance: ");
            scanf("%lf %lf %lf", &a, &b, &c);
            equationresult = bisection(a, b, c, functiondata);
            printf("Bisection's result is: %lf\n", equationresult);
            free(functiondata);
            break;
        case 2:
        	//DONE
            // Regula Falsi
            printf("\nRegula Falsi:\n");
            functiondata = functionallocator();
            printf("Enter lower bound, upper bound, tolerance, and maximum iterations: ");
            scanf("%lf %lf %lf %d", &a, &b, &c, &max);
            equationresult = regulaFalsi(a, b, functiondata, c, max);
            printf("Regula Falsi's result is: %lf\n", equationresult);
            free(functiondata);
            break;
        case 3:
            //DONE
            // Newton Raphson
            printf("\nNewton Raphson:\n");
            functiondata = functionallocator();
            printf("Enter initial guess, tolerance, and maximum iterations: ");
            scanf("%lf %lf %d", &a, &c, &max);
            equationresult = newtonRaphson(a, functiondata, c, max);
            printf("Newton Raphson's result is: %lf\n", equationresult);
            free(functiondata);
            break;
        case 4:
        	//DONE
            //Inverse Matrix
            printf("\nInverse Matrix:\n");
            inversematrixcall();
            break;
        case 5:
        	//DONE
            //Gauss Elimination
            printf("\nGauss Elimination:\n");
            gausseliminationcall();
            break;
        case 6:
        	//DONE
            //Gauss Seidal
            printf("\nGauss Seidal:\n");
            gaussSeidelCall();
            break;
        case 7:
        	//DONE
            // Numerical Derivative
            printf("\nNumerical Derivative:\n");
            numericalDerivativeCall();
            break;
        case 8:
        	//DONE
            // Simpson's 1/3 Method
            printf("\nSimpson's 1/3 Method:\n");
            simpsonOneThirdCall();
            break;
        case 9:
        	//DONE
            // Trapezoidal Method
            printf("\nTrapezoid Method:\n");
            trapezoidalRuleCall();
            break;
        case 10:
        	//DONE
            // Gregory Newton
            printf("\nGregory Newton:\n");
            gregoryNewtonBackwardCall();
            break;
    }
    return 0;
}
int main(){
    int userchoice;
    userchoice = menu();
    while (userchoice > 0 && userchoice < 11){
        methodchooser(userchoice);
        userchoice = menu();
    }
    return 0;
}

