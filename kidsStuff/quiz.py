# A multiple choice quiz

print("Hello, what is your name? ")
name = raw_input()
print "hello ", name, " would you like a quiz? (y/n)"

if raw_input() == 'n':
    print "ok, maybe another time, bye bye"
    exit()

score = 0

print "ok, here we go. Your score is currently ", score, " ..."

print "What has 4 legs, soft fur, pointy ears and says meeow?"
print "1 a dog"
print "2 a cat"
print "3 a snake"
print "4 a spikder"

ans = int(raw_input())

if ans == 2:
    print "Yes, that's right, it's a cat! Well done"
    score = score + 1
    print "Your score is now ", score
else:
    print "No silly, it's a cat."
    print "Your score is still ", score





