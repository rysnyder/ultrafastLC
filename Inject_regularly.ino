int outPin = 4;                 // LED connected to digital pin 4
int inPin = 7;                  // pushbutton connected to digital pin 7
int posAPin = 10;               // command position A (VICI P5)
int posBPin = 11;               // command position B (VICI P6)
int prev_val = 0;               // variable to store the read value

void setup() {
  pinMode(outPin, OUTPUT);      // LED output
  pinMode(posAPin, OUTPUT);     // Position A command output
  pinMode(posBPin, OUTPUT);     // Position B command output
  pinMode(inPin, INPUT_PULLUP); // pushbutton input
  Serial.begin(9600);           // initialize serial communication
  digitalWrite(outPin,HIGH);    // LED off at start
}

void loop() {
  digitalWrite(posAPin, LOW);   // command Position A 
  digitalWrite(posBPin, HIGH);  // ensure Position B not commanded
  digitalWrite(outPin, LOW);    // turn LED on during loading 
  delay(1500);                  // inject 1500 ms delay (for loading sample)
  digitalWrite(outPin, HIGH);   // turn LED off 
  digitalWrite(posAPin, HIGH);  // ensure Position A not commanded
  digitalWrite(posBPin, LOW);   // command Position B
  delay(4500);                  // injection frequency (4500 ms)
 
}                               // loop repeats



