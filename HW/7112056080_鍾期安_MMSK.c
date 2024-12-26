#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double arrivalTime;
    double serviceInterval;
    double serviceTime;
} Event;

typedef struct {
    Event *events;
    int capacity;
    int size;
    int front;
    int rear;
} Queue;

Queue* initializeQueue(int capacity) {
    Queue *queue = (Queue*)malloc(sizeof(Queue));
    queue->events = (Event*)malloc(sizeof(Event) * capacity);
    queue->capacity = capacity;
    queue->size = 0;
    queue->front = 0;
    queue->rear = -1;
    return queue;
}

void enqueue(Queue *queue, Event event) {
    if (queue->size < queue->capacity) {
        queue->rear = (queue->rear + 1) % queue->capacity;
        queue->events[queue->rear] = event;
        queue->size++;
    } else {
        printf("Queue is full.\n");
    }
}

Event dequeue(Queue *queue) {
    if (queue->size > 0) {
        Event event = queue->events[queue->front];
        queue->front = (queue->front + 1) % queue->capacity;
        queue->size--;
        return event;
    } else {
        //printf("Queue is empty.\n");
        Event event = {0.0, 0.0, 0.0};  // 空事件
        return event;
    }
}

int getQueueSize(Queue *queue) {
    return queue->size;
}

// 定義節點結構
typedef struct Node {
    Event event;
    struct Node* next;
} Node;

// 定義 linked list 結構
typedef struct {
    Node* head;
    size_t size;
    size_t maxSize;
} LinkedList;

// 初始化 linked list
void initializeLinkedList(LinkedList* list, size_t maxSize) {
    list->head = NULL;
    list->size = 0;
    list->maxSize = maxSize;
}

// 在 linked list 中插入事件，按照 serviceTime + serviceInterval 升冪排序
void insertEvent(LinkedList* list, Event newEvent) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    if (newNode == NULL) {
        printf("Memory allocation error.\n");
    }

    newNode->event = newEvent;
    newNode->next = NULL;

    // 插入新節點
    if (list->head == NULL || (newEvent.serviceTime + newEvent.serviceInterval) < (list->head->event.serviceTime + list->head->event.serviceInterval) ) {
        newNode->next = list->head;
        list->head = newNode;
    } else {
        Node* current = list->head;
        while (current->next != NULL && (newEvent.serviceTime + newEvent.serviceInterval) > (list->head->event.serviceTime + list->head->event.serviceInterval)) {
            current = current->next;
        }
        newNode->next = current->next;
        current->next = newNode;
    }

    // 更新 size
    list->size++;
}

// 釋放 linked list 佔用的記憶體
void freeLinkedList(LinkedList* list) {
    Node* current = list->head;
    while (current != NULL) {
        Node* temp = current;
        current = current->next;
        free(temp);
    }
}

// 回傳 linked list 中 DepartureTime 最小的值
double getMinDepartureTime(LinkedList* list) {
    if (list->head == NULL) {
        //printf("Linked list is empty.\n");
        return -1;
    }

    return list->head->event.serviceTime + list->head->event.serviceInterval;
}

// 刪除 linked list 中 serviceTime + serviceInterval 最小的節點
Event deleteMinDepartureTime(LinkedList* list) {

    Node* temp = list->head;
    list->head = list->head->next;
    Event deletedEvent = temp->event;
    free(temp);
    list->size--;

    return deletedEvent;
}

int factorial(int n) {
    if (n < 0) {
        return 0;
    } else if (n == 0 || n == 1) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

double generateExponentialRandom(double rate) {
    double u = rand() / (RAND_MAX + 1.0);
    return -log(1 - u) / rate;
}

void simulatedQueueSystem(double lambda, double mu, int s, int k) {
    double simulationTime = 100000.0;
    double clock = 0.0;
    Queue *queue = initializeQueue(k);
    double totalCustomersInQueue = 0.0;
    double totalCustomers = 0.0;
    double totalWaitingTime = 0.0;
    double totalWaitingTimeInQueue = 0.0;
    Event nextEvent, departureEvent, dequeueEvent;
    int customers = 0;
    LinkedList serviceDesk;
    initializeLinkedList(&serviceDesk, s);
    double arrivalInterval, serviceInterval;
    double nextArrivalTime;

    while (clock < simulationTime){
        
        arrivalInterval = generateExponentialRandom(lambda);
        serviceInterval = generateExponentialRandom(mu);
        nextEvent.arrivalTime = clock + arrivalInterval;
        nextEvent.serviceInterval = serviceInterval;
        nextArrivalTime = clock + arrivalInterval;

        while ( customers > 0 && getMinDepartureTime(&serviceDesk) <= nextArrivalTime ) {

            departureEvent = deleteMinDepartureTime(&serviceDesk);
            if ( getQueueSize(queue) > 0 ) {
                dequeueEvent = dequeue(queue);
                dequeueEvent.serviceTime = departureEvent.serviceTime + departureEvent.serviceInterval;
                insertEvent(&serviceDesk, dequeueEvent);
                totalCustomersInQueue += ((departureEvent.serviceTime + departureEvent.serviceInterval - clock) * (customers - s));
                totalWaitingTimeInQueue += (departureEvent.serviceTime + departureEvent.serviceInterval - departureEvent.arrivalTime);
            }
            totalCustomers += ((departureEvent.serviceTime + departureEvent.serviceInterval - clock) * customers);
            totalWaitingTime += (departureEvent.serviceTime + departureEvent.serviceInterval - departureEvent.arrivalTime);
            customers --;
            clock = departureEvent.serviceTime + departureEvent.serviceInterval;

        }
        
        totalCustomers += (nextArrivalTime - clock) * customers;

        if ( 0 <= customers < s ) {
            
            nextEvent.serviceTime = nextArrivalTime;
            insertEvent(&serviceDesk, nextEvent);
            customers ++;

        } else if ( s <= customers < (s + k) ) {

            nextEvent.serviceTime = 1000000;
            enqueue(queue, nextEvent);
            totalCustomersInQueue += (nextArrivalTime - clock) * (customers - s);
            customers ++;
            
        } else {
            totalCustomersInQueue += (nextArrivalTime - clock) * k;
        }

        clock = nextArrivalTime;
    }

    while ( customers > 0 && getMinDepartureTime(&serviceDesk) <= nextArrivalTime ) {

        departureEvent = deleteMinDepartureTime(&serviceDesk);
        if ( getQueueSize(queue) > 0 ) {
            dequeueEvent = dequeue(queue);
            dequeueEvent.serviceTime = departureEvent.serviceTime + departureEvent.serviceInterval;
            insertEvent(&serviceDesk, dequeueEvent);
            totalCustomersInQueue += ((departureEvent.serviceTime + departureEvent.serviceInterval - clock) * (customers - s));
            totalWaitingTimeInQueue += (departureEvent.serviceTime + departureEvent.serviceInterval - departureEvent.arrivalTime);
        }
        totalCustomers += ((departureEvent.serviceTime + departureEvent.serviceInterval - clock) * customers);
        totalWaitingTime += (departureEvent.serviceTime + departureEvent.serviceInterval - departureEvent.arrivalTime);
        customers --;
        clock = departureEvent.serviceTime + departureEvent.serviceInterval;

    }

    double L = totalCustomers / nextArrivalTime;  // 平均顧客數
    double Lq = totalCustomersInQueue / nextArrivalTime; // 平均隊列長度
    double W = totalWaitingTime / nextArrivalTime;  // 平均逗留时间
    double Wq = totalWaitingTimeInQueue / nextArrivalTime; // 平均等待时间

    printf("Simulated analysis:\n");
    printf("Lambda: %.2f, mu: %.1f, S: %d, K: %d\n", lambda, mu, s, k);
    printf("Average Number of Customers (L): %.2f\n", L);
    printf("Average Queue Length (Lq): %.3f\n", Lq);
    printf("Average Time in System (W): %.5f\n", W);
    printf("Average Time in Queue (Wq): %.5f\n", Wq);

    free(queue->events);
    free(queue);
    freeLinkedList(&serviceDesk);
}

double p_n(double lambda, double mu, int s, int k, int n){
    double Rho = lambda / ( s * mu );
    double p_0 = 1;
    for (int i = 1; i < s; i++) {
        p_0 += pow( ( lambda / mu ), i ) / factorial(i);
    }
    p_0 += (pow( ( lambda / mu ), s ) / factorial(s)) * ((1 - pow( Rho, k - s + 1 )) / ( 1 - Rho ));
    p_0 = 1 / p_0;

    if ( n == 0 ) {
        return p_0;
    }
    else if ( 0 < n < s ) {
        return (pow( ( lambda / mu ), n ) / factorial(n)) * p_0;
    } 
    else if ( s <= n <= k ) {
        return (pow( ( lambda / mu ), n ) / (factorial(s) * pow( s, n - s))) * p_0;
    } else {
        return 0;
    }
}

void theoreticalQueueSystem(double lambda, double mu, int s, int k) {

    double L = 0;
    double Lq = 0;
    double W = 0;
    double Wq = 0;
    double Rho = lambda / ( s * mu );
    double p_0 = p_n(lambda, mu, s, k, 0);
    double lambda_eff = lambda * ( 1 -  p_n(lambda, mu, s, k, k));

    Lq = (pow(( lambda / mu ), s ) * Rho ) / ( factorial(s) * ( 1 -  Rho ) * ( 1 -  Rho ));
    Lq *= p_0;
    Lq *= (( k - s ) * ( Rho - 1 ) * pow( Rho, k - s )) - pow( Rho, k - s ) + 1;

    L += Lq;
    for (int n = 0; n < s; n++) {
        L += n * p_n(lambda, mu, s, k, n);
    }
    L += s;
    for (int n = 0; n < s; n++) {
        L -= p_n(lambda, mu, s, k, n) * s;
    }

    W = L / lambda_eff;

    Wq = Lq / lambda_eff;
    
    printf("Theoretical analysis:\n");
    printf("Lambda: %.2f, mu: %.1f, S: %d, K: %d\n", lambda, mu, s, k);
    printf("Average Number of Customers (L): %.2f\n", L);
    printf("Average Queue Length (Lq): %.3f\n", Lq);
    printf("Average Time in System (W): %.5f\n", W);
    printf("Average Time in Queue (Wq): %.5f\n", Wq);

}

int main() {
    int s = 3;
    int k = 8;
    double mu = 40;

    double lambda;
    for (lambda = 15.0; lambda <= 25.0; lambda += 1.0) {
        theoreticalQueueSystem(lambda, mu, s, k);
        printf("\n");
        simulatedQueueSystem(lambda, mu, s, k);
        printf("\n");
        printf("////////////////////////////////////////////////////////////////\n");
    }

    return 0;
}