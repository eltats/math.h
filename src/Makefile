CC= gcc -std=c11 #-Wall -Wextra -Werror 
LIBSOURCES= $(wildcard s21_*.c)
LIBOBJECTS= $(LIBSOURCES:.c=.o)
COVERAGE= --coverage


all: test
	make -j gcov_report

test: s21_math.a
	$(CC) --coverage tests/tests.c -L. s21_math.a -lcheck `pkg-config --cflags --libs check`
	exec ./a.out

s21_math.a: $(LIBOBJECTS)
	ar rc s21_math.a $^
	ranlib s21_math.a

clean:
	rm -rf s21_math.a
	rm -rf *.o
	rm -rf a.out
	find . -name '*.gcno' -type f -delete
	find . -name '*.gcda' -type f -delete
	find . -name '*.info' -type f -delete
	rm -rf report

re: clean all

style:
	python3 ../materials/linters/cpplint.py --extensions=c,h *.c *.h

git:
	make clean
	git add .
	git commit -m 'test'
	git push

gcov_report: s21_math.a test
	lcov -t "test" -o test.info -c -d . --rc lcov_branch_coverage=1
	genhtml -o report test.info --rc lcov_branch_coverage=1

.PHONY: all clean test s21_math.a gcov_report
