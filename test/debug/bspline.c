#include <stdio.h>

int main() {
    /** Step 0 */
    int n, i, j;
    float xi,yi;
    printf("number of points:\n");
    scanf("%d", &n);
    
    n--;
    float x[n + 1], a[n + 1], h[n], A[n], l[n + 1],
        u[n + 1], z[n + 1], c[n + 1], b[n], d[n];
     printf("\n inter values for x\n");
    for (i = 0; i < n + 1; ++i) scanf("%f", &x[i]);
     printf("\n inter values for y\n");
    for (i = 0; i < n + 1; ++i) scanf("%f", &a[i]);

    /** Step 1 */
    for (i = 0; i <= n - 1; ++i) h[i] = x[i + 1] - x[i];

    /** Step 2 */
    for (i = 1; i <= n - 1; ++i)
        A[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

    /** Step 3 */
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    /** Step 4 */
    for (i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    /** Step 5 */
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    /** Step 6 */
    for (j = n - 1; j >= 0; --j) {
        c[j] = z[j] - u[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    /** Step 7 */
    printf("%2s %8s %8s %8s %8s\n", "i", "ai", "bi", "ci", "di");
    for (i = 0; i < n; ++i)
        printf("%2d %8.2f %8.2f %8.2f %8.2f\n", i, a[i], b[i], c[i], d[i]);
    printf("\n value needs to be interpolated \n");
    scanf("%lf",&xi);
    
    return 0;
}
