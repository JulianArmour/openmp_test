//
// Created by Julian on 2021-03-22.
//
#include <stdio.h>


int main(void) {
  struct A {int a;} a, *pa;
  a.a = 1;
  pa = &a;
  pa->a++;
  printf("%d", a.a);
}