def sendreport(mails = ["claudio.duran@ua.pt"],subject = "MRR report",content = "Daily report from MRR"):

    import smtplib, ssl

    class Mail:

        def __init__(self):
            self.port = 465
            self.smtp_server_domain_name = "smtp.gmail.com"
            self.sender_mail = "mrr.sender@gmail.com"
            self.password = "7f4nCx2fGqCSvwp"

        def send(self, emails, subject, content):
            ssl_context = ssl.create_default_context()
            service = smtplib.SMTP_SSL(self.smtp_server_domain_name, self.port, context=ssl_context)
            service.login(self.sender_mail, self.password)
            



            for email in emails:
                result = service.sendmail(self.sender_mail, email, f"Subject: {subject}\n{content}")

            service.quit()



    mail = Mail()

    print("Sending report to :", mails)

    mail.send(mails, subject, content)



