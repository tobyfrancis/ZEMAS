#include <SoftwareSerial.h>

SoftwareSerial Serial3(9,10); // RX, TX

char buffer[10] = {0};
char pos[10] = {0};
char inc[10] = {0};
char text[10] = {0};
char text2[10] = {0};
char text3[10] = {0};
char runTo[64];
char command[64];

// Define PINS
int coarse = 2;
int fine = 3;
int up = 4;
int down = 5;
// int Analog_Out = 14;
int Analog_Signal = 15;
int joystick_pin = 14;
int an_val = 0;

// Define variables
float currentValue = 0;
float increment = 10;
String textfield;
String str1 = String("step_size");
String str2 = String("stop_interval");
String str3 = String("RPM");
String str4 = String("Duration");
String str5 = String("Iterations");
int enterState = 0;
float val = 0.0;
int num = 0;
int dotState = 0;
int decimal = 0;
int mode  =  0;
float step_size = 0;
float stop_interval = 0;
int RPM = 0;
float resolution = 3420160;
long translated_angle = 0;
float calculated_position = 0;
int negative = 0;
int dec_serial = 0;
int char_serial = 0;
int calculate = 1;
int raw_position = 0;
float position_old = 0;
int first = 0;
int iterations_val = 0;
int bigNumber = 0;
int rotation = 0;
float mult_inc = 0;
int bias_state = 0;


// Funtions for controller
void Initialize_Motor()
{
  Serial3.println("EN");
  Serial3.println("SOR0");
  Serial3.println("CONTMOD");
  // Serial3.println("SP4000");
}

void Start_Motor()
{
  Serial3.println("SR120");
  Serial3.println("NP");  // Notify Position
  Serial3.println("M"); // Initiate motion
  delay(5);
  while (Serial3.read()!=112) // 112 is Dec for "p"
  {
    // Serial.println("STUCK");
    delay(1);
  }
  Serial.print(currentValue);
  Serial.print(": DONE");
  Serial.print("\n");
}

float Update_Position(int first)
{
  calculated_position = 0;
  negative = 0;
  dec_serial = 0;
  char_serial = 0;
  Serial3.println("POS");
  delay(15);
  while (Serial3.available())
  {
    calculate = 1;
    dec_serial = Serial3.read();
    char_serial = dec_serial - 48;
    // If Line feed and Carriage Feed
    if (dec_serial == 13 || dec_serial == 10)
    {
      calculate = 0;
      char_serial = 0;
    }
    // If the char number is not between 0 and 9 then garbage
    // no need to calculate
    if (dec_serial <45 || dec_serial>57)
    {
      char_serial = 0;
      calculate = 0;
    }

    // If negative value is returned
    if (dec_serial == 45)
    {
      negative = 1;
      char_serial = 0;
      calculate = 0;
    }

    // If positive, not garbage then we can use the data
    if (calculate==1)
    {
      raw_position = raw_position + char_serial;
      raw_position = raw_position*10;  // why *10?
    }
    delay(1);
  }
  raw_position = raw_position/10;
  if (negative == 1)
  {
    raw_position = raw_position*-1;
  }
  calculated_position = (float)raw_position/3420160*360;
  if (first==1)
  {
    calculated_position==0;
  }
  return calculated_position;
}



// Page 1
void biasPopCallback(void *ptr)
{
    page2.show();
    bias_state = 1;
}

void tomoPopCallback(void *ptr)
{
    page3.show();
    bias_state = 0;
}


// Page 2
void fine_buttonPopCallback(void *ptr)
{
    increment -= 0.01;
    dtostrf(increment, 8, 2, inc);
    t3.setText(inc);
}

void coarse_buttonPopCallback(void *ptr)
{
    increment += 0.01;
    dtostrf(increment, 8, 2, inc);
    t3.setText(inc);
}

void up_buttonPopCallback(void *ptr)
{
    currentValue = currentValue + increment;
    dtostrf(currentValue, 7, 2, pos);
    t2.setText(pos);
    // Sending command
    // Initialize_Motor();
    translated_angle = (increment*resolution/360.0);
    sprintf(command, "LR%ld", translated_angle);  // Load Relative Position
    Serial.print(command);
    Serial.print("\n");
    Serial3.println(command);
    Start_Motor();
    calculated_position = Update_Position(0);
    // dtostrf(calculated_position, 7,3,pos);
    // t2.setText(pos);
    position_old = calculated_position;
}

void down_buttonPopCallback(void *ptr)
{
  currentValue = currentValue - increment;
  dtostrf(currentValue, 7, 2, pos);
  t2.setText(pos);
  // Sending command
  Initialize_Motor();
  translated_angle = (increment*resolution/360);
  sprintf(command, "LR-%ld", translated_angle);
  Serial.print(command);
  Serial.print("\n");
  Serial3.println(command);
  Start_Motor();
  calculated_position = Update_Position(0);
  position_old = calculated_position;
}


void back_page2PopCallback(void *ptr)
{
  page1.show();
  bias_state = 0;
}


// Page 3
void mode1PopCallback(void *ptr)
{
  bigNumber = 0;
  page4.show();
  Serial3.println("HO");
  mode = 1;
}

void mode2PopCallback(void *ptr)
{
  bigNumber = 1;
  page5.show();
  mode = 2;
}

void mode3PopCallback(void *ptr)
{
  bigNumber = 1;
  page7.show();
  mode = 3;
}

void back_page3PopCallback(void *ptr)
{
  page1.show();
}

// Page 4
void back_page4PopCallback(void *ptr)
{
  page3.show();
}

void step_sizeButtonPopCallback(void *ptr)
{
  textfield = String("step_size");
  val = 0;
  decimal = 0;
  dotState = 0;
  page6.show();
}

void stop_intervalButtonPopCallback(void *ptr)
{
  textfield = String("stop_interval");
  val = 0;
  decimal = 0;
  dotState = 0;
  page6.show();
}

void run_page4PopCallback(void *ptr)
{
  // Sending command to motor
  Serial3.println("EN");
  for (int i = 1; i<=iterations_val; i++)
  {
    translated_angle = step_size*resolution/360.0;
    sprintf(command, "LR%6ld", translated_angle);
    Serial.println(command);
    Serial3.println(command);
    Serial3.println("NP");
    Serial3.println("M");

    while (Serial3.read()!=112)
    {
      delay(1);
    }
    Serial.print("ONE ITERATION DONE\n");
    delay(stop_interval*1000);
  }
  Serial.println("OUT OF LOOP\n");

}

void iterationsPopCallback(void *ptr)
{
  textfield = String("Iterations");
  val = 0;
  page6.show();
}


// Page 5
void back_page5PopCallback(void *ptr)
{
  page3.show();
}

void RPM_buttonPopCallback(void *ptr)
{
  textfield = String("RPM");
  val = 0;
  decimal = 0;
  dotState = 0;
  // bigNumber = 1;
  page6.show();
}

void run_page5PopCallback(void *ptr)
{
  // Sending command to motor
  Serial3.println("NV");
  translated_angle = RPM*1670;
  sprintf(command, "v%d", translated_angle);
  Serial.println(command);
  Serial3.println(command); 
}

void stop_page5PopCallback(void *ptr)
{
  Serial3.println("v0");
  Serial.println("STOP\n");
}

// Page 6
void zeroPopCallback(void *ptr)
{
  num = 0;
  Display(textfield, num); 
}

void onePopCallback(void *ptr)
{
  num = 1;
  Display(textfield, num);
}

void twoPopCallback(void *ptr)
{
  num = 2;
  Display(textfield, num);
}


void threePopCallback(void *ptr)
{
  num = 3;
  Display(textfield, num);
}

void fourPopCallback(void *ptr)
{
  num = 4;
  Display(textfield, num);
}

void fivePopCallback(void *ptr)
{
  num = 5;
  Display(textfield, num);
}

void sixPopCallback(void *ptr)
{
  num = 6;
  Display(textfield, num);
}

void sevenPopCallback(void *ptr)
{
  num = 7;
  Display(textfield, num);
}

void eightPopCallback(void *ptr)
{
  num = 8;
  Display(textfield, num);
}

void ninePopCallback(void *ptr)
{
  num = 9;
  Display(textfield, num);
}

void Display(String textfield, int num)
{
  display_input.setText(KeyboardLogic1(num));
}


void enterPopCallback(void *ptr)
{
  enterState = 1;
  switch (mode)
  {
    case 1:
      page4.show();
      if (textfield==str1)
      {
        step_size = val;
        dtostrf(step_size, 8,2, text);
        dtostrf(stop_interval, 8,2, text2);
        sprintf(text3, "%d", iterations_val);
        display_angle.setText(text);
        display_interval.setText(text2);
        display_iterations.setText(text3);
      }
      else if (textfield==str2)
      {
        stop_interval = val;
        dtostrf(step_size, 8,2, text);
        dtostrf(stop_interval, 8,2, text2);
        sprintf(text3, "%d", iterations_val);
        display_angle.setText(text);
        display_interval.setText(text2);
        display_iterations.setText(text3);
      }
      else if (textfield==str5)
      {
        iterations_val = val;
        dtostrf(step_size, 8,2, text);
        dtostrf(stop_interval, 8,2, text2);
        sprintf(text3, "%d", iterations_val);
        display_angle.setText(text);
        display_interval.setText(text2);
        display_iterations.setText(text3);
      }
      break;
    case 2:
      page5.show();
      if (textfield==str3)
      {
        RPM = (int) val;
        Serial.println(RPM);
        sprintf(text, "%d", RPM);
        Serial.println(text);
        display_RPM.setText(text);
      }
      break;
    case 3:
      page7.show();
      rotation = (int) val;
      Serial.println(rotation);
      sprintf(text, "%d", rotation);
      Serial.println(text);
      display_rotation.setText(text);
      break;
  }

}


void delete_buttonPopCallback(void *ptr)
{
  val = 0;
  decimal = 0;
  dtostrf(val, 8, 2, text);
  display_input.setText(text);
}

void dotPopCallback(void *ptr)
{
  dotState = !dotState;
}

char* KeyboardLogic1(int num)
{   
  switch (bigNumber)
  {
  case 0:
    if (dotState==0)
    {
      if (val<10)
      {
        val = num + val*10;
      }
      else if (val==0)
      {
        val = num + val;
      }
      dtostrf(val, 8,0, text);
    }
    else if (dotState ==1)
    {
      if (decimal == 0)
      {
        val = val + (float)num/10;
        dtostrf(val, 8,1, text);
        decimal++;
      }
      else if (decimal == 1)
      {
        val = val + (float)num/100;
        dtostrf(val, 8,2, text);
        decimal++;
      } 
    }
    dtostrf(val, 8,2, text);
    break;

    case 1:
      val = val+num;
      int val2 = (int) val;
      sprintf(text, "%d", val2);
      Serial.println(val2);
      Serial.println(text);
      break;
    }
    Serial.print(text);
    Serial.print("\n");  
    return text;
 
}

// Page 7
void rotation_buttonPopCallback(void *ptr)
{
  textfield = String("Rotation");
  val = 0;
  decimal = 0;
  dotState = 0;
  page6.show();
}

void back_page7PopCallback(void *ptr)
{
  page3.show();
}

void run_page7PopCallback(void *ptr)
{
  translated_angle = rotation*3420160;
  sprintf(command, "%d", translated_angle);
  Serial3.println("NP");
  Serial3.println(command);
  Serial3.println("M");
}




void setup() {
  nexInit();
  // Page 1
  bias.attachPop(biasPopCallback, &bias);
  tomo.attachPop(tomoPopCallback, &tomo);
  // Page 2
  fine_button.attachPop(fine_buttonPopCallback, &fine_button);
  coarse_button.attachPop(coarse_buttonPopCallback, &coarse_button);
  up_button.attachPop(up_buttonPopCallback, &up_button);
  down_button.attachPop(down_buttonPopCallback, &down_button);
  back_page2.attachPop(back_page2PopCallback, &back_page2);
  // Page 3
  mode1.attachPop(mode1PopCallback, &mode1);
  mode2.attachPop(mode2PopCallback, &mode2);
  mode3.attachPop(mode3PopCallback, &mode3);
  back_page3.attachPop(back_page3PopCallback, &back_page3);
  // Page 4
  step_sizeButton.attachPop(step_sizeButtonPopCallback, &step_sizeButton);
  stop_intervalButton.attachPop(stop_intervalButtonPopCallback, &stop_intervalButton);
  back_page4.attachPop(back_page4PopCallback, &back_page4);
  run_page4.attachPop(run_page4PopCallback, &run_page4);
  iterations.attachPop(iterationsPopCallback, &iterations);
  // Page 5
  RPM_button.attachPop(RPM_buttonPopCallback, &RPM_button);
  stop_page5.attachPop(stop_page5PopCallback, &stop_page5);
  back_page5.attachPop(back_page5PopCallback, &back_page5);
  run_page5.attachPop(run_page5PopCallback, &run_page5);
  // Page 6 
  zero.attachPop(zeroPopCallback, &zero);
  one.attachPop(onePopCallback, &one);
  two.attachPop(twoPopCallback, &two);
  three.attachPop(threePopCallback, &one);
  four.attachPop(fourPopCallback, &one);
  five.attachPop(fivePopCallback, &one);
  six.attachPop(sixPopCallback, &one);
  seven.attachPop(sevenPopCallback, &one);
  eight.attachPop(eightPopCallback, &one);
  nine.attachPop(ninePopCallback, &one);
  delete_button.attachPop(delete_buttonPopCallback, &delete_button);
  dot.attachPop(dotPopCallback, &dot);
  enter.attachPop(enterPopCallback, &enter);
  // Page 7
  rotation_button.attachPop(rotation_buttonPopCallback, &rotation_button);
  back_page7.attachPop(back_page7PopCallback, &back_page7);
  run_page7.attachPop(run_page7PopCallback, &run_page7);
  
  Serial.begin(9600);

  pinMode(coarse, INPUT_PULLUP);
  pinMode(fine, INPUT_PULLUP);
  pinMode(up, INPUT_PULLUP);
  pinMode(down, INPUT_PULLUP);

  // pinMode(Analog_Out, OUTPUT);
  pinMode(Analog_Signal, OUTPUT);
  // digitalWrite(Analog_Out, LOW);
  digitalWrite(Analog_Signal, LOW);
  pinMode(joystick_pin, INPUT);


  Serial3.begin(9600);
  Serial3.println("EN");  // Enable Drive
  Serial3.println("SOR0");  //Source for velocity (0 for interface)
  Serial3.println("CONTMOD"); // Continuous mode
  Serial3.println("MAV50");
  Serial3.println("SP5000");
  Serial3.println("PP63");  // Load Position Proportional Term
  Serial3.println("I2");  
  Serial3.println("POR12");
  Serial3.println("PD1"); // Load Position Differential Term
  Serial3.println("CORRIDOR40");
  Serial3.println("HO");  // Define Home Position
  Serial3.println("LL-999999999");  // Setting Lower Limit
  Serial3.println("LL999999999");   // Setting Upper Limit
  Serial.print("Setup is complete\n");
}

void loop() {
  nexLoop(nex_Listen_List);
  if (bias_state ==1)
  {
    an_val = analogRead(joystick_pin);
    // Serial.println(an_val);
    if (an_val<=420)
    {
      Serial.println("LEFT x 1");
      // increment: 0.01
      mult_inc = increment*1;
      currentValue = currentValue + mult_inc;
      dtostrf(currentValue, 7,2, pos);
      t2.setText(pos); 
      translated_angle = (mult_inc*resolution/360.0);
      sprintf(command, "LR%ld", translated_angle);
      Serial.println(command);
      Serial3.println(command);
      Start_Motor();
    }

    if (an_val>=560)
    {
      Serial.println("RIGHT x 1");
      // increment: 0.01
      mult_inc = increment*1;
      currentValue = currentValue - mult_inc;
      dtostrf(currentValue, 7,2, pos);
      t2.setText(pos); 
      translated_angle = (mult_inc*resolution/360.0);
      sprintf(command, "LR-%ld", translated_angle);
      Serial.println(command);
      Serial3.println(command);
      Start_Motor();
    }
  }

}
